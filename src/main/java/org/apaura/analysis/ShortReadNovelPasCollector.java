package org.apaura.analysis;

import org.apaura.io.AlignmentReader;
import org.apaura.io.AlignmentReaderFactory;
import org.apaura.io.AlignmentRecord;
import org.apaura.model.Interval;
import org.apaura.model.Sample;
import org.apaura.model.Strand;
import org.apaura.model.TranscriptModel;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/**
 * Collects CIGAR-aware terminal-exon coverage from conventional RNA-seq
 * alignments and performs condition-blind novel PAS breakpoint inference.
 *
 * Sample condition labels are deliberately never used here.
 */
public final class ShortReadNovelPasCollector {

    private static final int BIN_SIZE = 16_384;

    private ShortReadNovelPasCollector() {
    }

    public static final class Site {
        public final String familyKey;
        public final String geneId;
        public final String geneName;
        public final String chrom;
        public final Strand strand;
        public final int position0;
        public final int supportingSamples;
        public final double upstreamDepth;
        public final double downstreamDepth;
        public final double coverageDrop;
        public final int distanceToAnnotatedPas;
        public final String status;

        Site(
                String familyKey,
                String geneId,
                String geneName,
                String chrom,
                Strand strand,
                int position0,
                int supportingSamples,
                double upstreamDepth,
                double downstreamDepth,
                double coverageDrop,
                int distanceToAnnotatedPas,
                String status
        ) {
            this.familyKey = familyKey;
            this.geneId = geneId;
            this.geneName = geneName;
            this.chrom = chrom;
            this.strand = strand;
            this.position0 = position0;
            this.supportingSamples = supportingSamples;
            this.upstreamDepth = upstreamDepth;
            this.downstreamDepth = downstreamDepth;
            this.coverageDrop = coverageDrop;
            this.distanceToAnnotatedPas = distanceToAnnotatedPas;
            this.status = status;
        }

        public boolean accepted() {
            return "accepted".equals(status);
        }
    }

    public static final class Result {
        public final List<Site> sites;
        public final Map<String, List<Integer>> acceptedByFamily;

        Result(
                List<Site> sites,
                Map<String, List<Integer>> acceptedByFamily
        ) {
            this.sites = Collections.unmodifiableList(
                    new ArrayList<Site>(sites)
            );

            Map<String, List<Integer>> copied =
                    new LinkedHashMap<String, List<Integer>>();

            for (Map.Entry<String, List<Integer>> entry
                    : acceptedByFamily.entrySet()) {

                copied.put(
                        entry.getKey(),
                        Collections.unmodifiableList(
                                new ArrayList<Integer>(
                                        entry.getValue()
                                )
                        )
                );
            }

            this.acceptedByFamily =
                    Collections.unmodifiableMap(copied);
        }
    }

    private static final class Target {
        final TerminalFamilyReconciler.Family family;
        final String key;
        final int start0;
        final int end0;

        Target(
                TerminalFamilyReconciler.Family family,
                int start0,
                int end0
        ) {
            this.family = family;
            this.key = familyKey(family);
            this.start0 = start0;
            this.end0 = end0;
        }

        int length() {
            return end0 - start0;
        }
    }

    private static final class TargetIndex {
        final Map<String, Map<Integer, List<Target>>> bins =
                new HashMap<String, Map<Integer, List<Target>>>();

        TargetIndex(List<Target> targets) {
            for (Target target : targets) {
                int firstBin = target.start0 / BIN_SIZE;
                int lastBin = (target.end0 - 1) / BIN_SIZE;

                Map<Integer, List<Target>> chromosomeBins =
                        bins.computeIfAbsent(
                                target.family.chrom,
                                key -> new HashMap<Integer, List<Target>>()
                        );

                for (int bin = firstBin; bin <= lastBin; bin++) {
                    chromosomeBins.computeIfAbsent(
                            bin,
                            key -> new ArrayList<Target>()
                    ).add(target);
                }
            }
        }

        Set<Target> candidates(
                String chrom,
                int start0,
                int end0
        ) {
            Map<Integer, List<Target>> chromosomeBins =
                    bins.get(chrom);

            if (chromosomeBins == null) {
                return Collections.emptySet();
            }

            Set<Target> targets =
                    new HashSet<Target>();

            int firstBin = Math.max(
                    0,
                    start0 / BIN_SIZE
            );

            int lastBin = Math.max(
                    0,
                    (end0 - 1) / BIN_SIZE
            );

            for (int bin = firstBin; bin <= lastBin; bin++) {
                List<Target> values =
                        chromosomeBins.get(bin);

                if (values != null) {
                    targets.addAll(values);
                }
            }

            return targets;
        }
    }

    public static Result discover(
            List<Sample> samples,
            List<TranscriptModel> transcripts,
            int anchorTolerance,
            int pasMergeDistance,
            int minMapq,
            boolean keepDuplicates,
            String library,
            int window,
            double minDepth,
            double minDrop,
            int minReplicates,
            int clusterDistance,
            int minSegmentLength
    ) throws IOException {

        List<Target> targets = buildTargets(
                transcripts,
                anchorTolerance,
                pasMergeDistance,
                minSegmentLength
        );

        if (targets.isEmpty()) {
            return new Result(
                    Collections.<Site>emptyList(),
                    Collections.<String, List<Integer>>emptyMap()
            );
        }

        TargetIndex index = new TargetIndex(targets);

        Map<String, List<ShortReadPasDiscovery.SampleCall>>
                callsByFamily =
                new LinkedHashMap<
                        String,
                        List<ShortReadPasDiscovery.SampleCall>
                >();

        /*
         * Process one biological sample at a time. Genome-wide depth arrays
         * are discarded after breakpoint evidence is extracted.
         */
        for (Sample sample : samples) {

            Map<String, int[]> deltaByFamily =
                    new LinkedHashMap<String, int[]>();

            for (Target target : targets) {
                deltaByFamily.put(
                        target.key,
                        new int[target.length() + 1]
                );
            }

            try (AlignmentReader reader =
                         AlignmentReaderFactory.open(sample.path)) {

                AlignmentRecord alignment;

                while ((alignment = reader.next()) != null) {
                    if (!isUsable(
                            alignment,
                            minMapq,
                            keepDuplicates
                    )) {
                        continue;
                    }

                    Strand transcriptStrand =
                            transcriptStrand(
                                    alignment,
                                    library
                            );

                    for (Interval block : alignment.alignedBlocks) {

                        for (Target target : index.candidates(
                                alignment.referenceName,
                                block.start0,
                                block.end0
                        )) {
                            if (!strandMatches(
                                    target.family.strand,
                                    transcriptStrand,
                                    library
                            )) {
                                continue;
                            }

                            int overlapStart = Math.max(
                                    block.start0,
                                    target.start0
                            );

                            int overlapEnd = Math.min(
                                    block.end0,
                                    target.end0
                            );

                            if (overlapEnd <= overlapStart) {
                                continue;
                            }

                            int[] delta =
                                    deltaByFamily.get(target.key);

                            int localStart =
                                    overlapStart - target.start0;

                            int localEnd =
                                    overlapEnd - target.start0;

                            delta[localStart]++;
                            delta[localEnd]--;
                        }
                    }
                }
            }

            for (Target target : targets) {
                int[] delta =
                        deltaByFamily.get(target.key);

                int[] depth =
                        new int[target.length()];

                int running = 0;

                for (int i = 0; i < depth.length; i++) {
                    running += delta[i];
                    depth[i] = running;
                }

                ShortReadPasDiscovery.SampleProfile profile =
                        new ShortReadPasDiscovery.SampleProfile(
                                sample.name,
                                target.family.chrom,
                                target.family.strand,
                                target.start0,
                                target.end0,
                                depth
                        );

                List<ShortReadPasDiscovery.SampleCall> calls =
                        ShortReadPasDiscovery.callSample(
                                profile,
                                window,
                                minDepth,
                                minDrop,
                                clusterDistance,
                                minSegmentLength
                        );

                if (!calls.isEmpty()) {
                    callsByFamily.computeIfAbsent(
                            target.key,
                            key -> new ArrayList<
                                    ShortReadPasDiscovery.SampleCall
                                    >()
                    ).addAll(calls);
                }
            }
        }

        List<Site> sites =
                new ArrayList<Site>();

        Map<String, List<Integer>> accepted =
                new LinkedHashMap<String, List<Integer>>();

        for (Target target : targets) {

            List<ShortReadPasDiscovery.SampleCall> calls =
                    callsByFamily.get(target.key);

            if (calls == null || calls.isEmpty()) {
                continue;
            }

            ShortReadPasDiscovery.Result result =
                    ShortReadPasDiscovery.aggregateCalls(
                            calls,
                            minReplicates,
                            clusterDistance
                    );

            for (ShortReadPasDiscovery.Candidate candidate
                    : result.candidates) {

                int distance =
                        distanceToAnnotatedPas(
                                candidate.position0,
                                target.family.pasPositions
                        );

                String status = candidate.status;

                if ("accepted".equals(status)) {
                    if (distance <= pasMergeDistance) {
                        status = "near_annotated_pas";
                    } else if (distance < minSegmentLength) {
                        status = "too_close_to_annotated_pas";
                    }
                }

                Site site = new Site(
                        target.key,
                        target.family.geneId,
                        target.family.geneName,
                        target.family.chrom,
                        target.family.strand,
                        candidate.position0,
                        candidate.supportingSamples,
                        candidate.upstreamDepth,
                        candidate.downstreamDepth,
                        candidate.coverageDrop,
                        distance,
                        status
                );

                sites.add(site);

                if (site.accepted()) {
                    accepted.computeIfAbsent(
                            target.key,
                            key -> new ArrayList<Integer>()
                    ).add(site.position0);
                }
            }
        }

        Collections.sort(
                sites,
                new Comparator<Site>() {
                    @Override
                    public int compare(
                            Site left,
                            Site right
                    ) {
                        int comparison =
                                left.chrom.compareTo(right.chrom);

                        if (comparison != 0) {
                            return comparison;
                        }

                        comparison = Integer.compare(
                                left.position0,
                                right.position0
                        );

                        if (comparison != 0) {
                            return comparison;
                        }

                        return left.geneId.compareTo(right.geneId);
                    }
                }
        );

        for (List<Integer> positions : accepted.values()) {
            Collections.sort(positions);
        }

        return new Result(
                sites,
                accepted
        );
    }

    static String familyKey(
            TerminalFamilyReconciler.Family family
    ) {
        return family.geneId
                + "\u0001"
                + family.chrom
                + "\u0001"
                + family.strand.symbol()
                + "\u0001"
                + family.anchor0;
    }

    public static void writeDiagnostics(
            File file,
            Result result
    ) throws IOException {

        try (PrintWriter writer =
                     new PrintWriter(
                             new BufferedWriter(
                                     new FileWriter(file)
                             )
                     )) {

            writer.println(
                    "gene_id\tgene_name\tchrom\tstrand\t"
                            + "position_1based\tsite_origin\t"
                            + "discovery_method\tsupporting_samples\t"
                            + "upstream_depth\tdownstream_depth\t"
                            + "coverage_drop\t"
                            + "distance_to_annotated_pas\tstatus"
            );

            for (Site site : result.sites) {
                writer.printf(
                        Locale.ROOT,
                        "%s\t%s\t%s\t%c\t%d\t"
                                + "novel\tcoverage_breakpoint\t"
                                + "%d\t%.6f\t%.6f\t%.6f\t%d\t%s%n",
                        site.geneId,
                        site.geneName,
                        site.chrom,
                        site.strand.symbol(),
                        site.position0 + 1,
                        site.supportingSamples,
                        site.upstreamDepth,
                        site.downstreamDepth,
                        site.coverageDrop,
                        site.distanceToAnnotatedPas,
                        site.status
                );
            }
        }
    }

    private static List<Target> buildTargets(
            List<TranscriptModel> transcripts,
            int anchorTolerance,
            int pasMergeDistance,
            int minSegmentLength
    ) {
        Map<String, List<TranscriptModel>> byGene =
                new LinkedHashMap<
                        String,
                        List<TranscriptModel>
                >();

        for (TranscriptModel transcript : transcripts) {
            byGene.computeIfAbsent(
                    transcript.geneId,
                    key -> new ArrayList<TranscriptModel>()
            ).add(transcript);
        }

        List<Target> targets =
                new ArrayList<Target>();

        for (List<TranscriptModel> geneTranscripts
                : byGene.values()) {

            TerminalFamilyReconciler.Result reconciliation =
                    TerminalFamilyReconciler.reconcile(
                            geneTranscripts,
                            anchorTolerance,
                            pasMergeDistance
                    );

            for (TerminalFamilyReconciler.Family family
                    : reconciliation.families) {

                if (family.pasPositions.isEmpty()) {
                    continue;
                }

                int minPas = Integer.MAX_VALUE;
                int maxPas = Integer.MIN_VALUE;

                for (Integer position : family.pasPositions) {
                    minPas = Math.min(minPas, position);
                    maxPas = Math.max(maxPas, position);
                }

                int start0;
                int end0;

                if (family.strand == Strand.MINUS) {
                    start0 = minPas;
                    end0 = family.anchor0;
                } else {
                    start0 = family.anchor0;
                    end0 = maxPas + 1;
                }

                if (end0 <= start0) {
                    continue;
                }

                if (end0 - start0
                        < 2 * minSegmentLength + 1) {
                    continue;
                }

                targets.add(
                        new Target(
                                family,
                                start0,
                                end0
                        )
                );
            }
        }

        return targets;
    }

    private static int distanceToAnnotatedPas(
            int position0,
            List<Integer> annotated
    ) {
        int best = Integer.MAX_VALUE;

        for (Integer position : annotated) {
            best = Math.min(
                    best,
                    Math.abs(position0 - position)
            );
        }

        return best;
    }

    private static boolean isUsable(
            AlignmentRecord alignment,
            int minMapq,
            boolean keepDuplicates
    ) {
        if (alignment.isUnmapped()
                || alignment.isSecondary()
                || alignment.isSupplementary()
                || alignment.isQcFail()
                || alignment.mapq < minMapq) {
            return false;
        }

        if (!keepDuplicates
                && alignment.isDuplicate()) {
            return false;
        }

        return !alignment.alignedBlocks.isEmpty();
    }

    private static boolean strandMatches(
            Strand targetStrand,
            Strand transcriptStrand,
            String library
    ) {
        if ("unstranded".equalsIgnoreCase(library)) {
            return true;
        }

        return transcriptStrand == Strand.UNKNOWN
                || targetStrand == transcriptStrand;
    }

    private static Strand transcriptStrand(
            AlignmentRecord alignment,
            String mode
    ) {
        if (mode == null
                || "unstranded".equalsIgnoreCase(mode)) {
            return Strand.UNKNOWN;
        }

        Strand alignmentStrand =
                alignment.alignmentStrand();

        if (isForwardMode(mode)) {
            return alignment.isPaired()
                    && alignment.isSecondOfPair()
                    ? alignmentStrand.flip()
                    : alignmentStrand;
        }

        if (isReverseMode(mode)) {
            return alignment.isPaired()
                    && alignment.isSecondOfPair()
                    ? alignmentStrand
                    : alignmentStrand.flip();
        }

        throw new IllegalArgumentException(
                "Unknown --library mode: "
                        + mode
                        + " (use unstranded, forward, or reverse)"
        );
    }

    private static boolean isForwardMode(
            String mode
    ) {
        return "forward".equalsIgnoreCase(mode)
                || "secondstrand".equalsIgnoreCase(mode)
                || "fr-secondstrand".equalsIgnoreCase(mode);
    }

    private static boolean isReverseMode(
            String mode
    ) {
        return "reverse".equalsIgnoreCase(mode)
                || "firststrand".equalsIgnoreCase(mode)
                || "fr-firststrand".equalsIgnoreCase(mode);
    }
}
