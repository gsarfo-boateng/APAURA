package org.apaura.analysis;

import org.apaura.io.AlignmentReader;
import org.apaura.io.AlignmentReaderFactory;
import org.apaura.io.AlignmentRecord;
import org.apaura.model.CoverageSegment;
import org.apaura.model.Interval;
import org.apaura.model.ShortReadGroup;
import org.apaura.model.ShortReadQuantification;
import org.apaura.model.Strand;
import org.apaura.stats.IsotonicRegression;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/** Quantifies annotated tandem terminal-exon APA from short-read RNA-seq coverage. */
public final class ShortReadQuantifier {
    private static final int BIN_SIZE = 16_384;

    private ShortReadQuantifier() {
        // Utility class.
    }

    private static final class SegmentRef {
        final int groupIndex;
        final int segmentIndex;
        final int start0;
        final int end0;
        final ShortReadGroup group;

        SegmentRef(
                int groupIndex,
                int segmentIndex,
                ShortReadGroup group,
                CoverageSegment segment
        ) {
            this.groupIndex = groupIndex;
            this.segmentIndex = segmentIndex;
            this.group = group;
            this.start0 = segment.start0;
            this.end0 = segment.end0;
        }
    }

    /** Simple genomic bin index used to avoid scanning every APA segment per read. */
    private static final class SegmentIndex {
        private final Map<String, Map<Integer, List<SegmentRef>>> bins =
                new HashMap<String, Map<Integer, List<SegmentRef>>>();

        SegmentIndex(List<ShortReadGroup> groups) {
            for (int groupIndex = 0; groupIndex < groups.size(); groupIndex++) {
                ShortReadGroup group = groups.get(groupIndex);
                addGroup(groupIndex, group);
            }
        }

        Set<SegmentRef> candidates(String chromosome, int start0, int end0) {
            Map<Integer, List<SegmentRef>> chromosomeBins = bins.get(chromosome);
            if (chromosomeBins == null) {
                return Collections.emptySet();
            }

            Set<SegmentRef> candidates = new HashSet<SegmentRef>();
            int firstBin = Math.max(0, start0 / BIN_SIZE);
            int lastBin = Math.max(0, (end0 - 1) / BIN_SIZE);

            for (int bin = firstBin; bin <= lastBin; bin++) {
                List<SegmentRef> refs = chromosomeBins.get(bin);
                if (refs != null) {
                    candidates.addAll(refs);
                }
            }
            return candidates;
        }

        private void addGroup(int groupIndex, ShortReadGroup group) {
            for (int segmentIndex = 0;
                    segmentIndex < group.segments.size();
                    segmentIndex++) {
                CoverageSegment segment = group.segments.get(segmentIndex);
                SegmentRef ref = new SegmentRef(
                        groupIndex,
                        segmentIndex,
                        group,
                        segment
                );

                int firstBin = segment.start0 / BIN_SIZE;
                int lastBin = (segment.end0 - 1) / BIN_SIZE;
                Map<Integer, List<SegmentRef>> chromosomeBins =
                        bins.computeIfAbsent(
                                group.chrom,
                                key -> new HashMap<Integer, List<SegmentRef>>()
                        );

                for (int bin = firstBin; bin <= lastBin; bin++) {
                    chromosomeBins.computeIfAbsent(
                            bin,
                            key -> new ArrayList<SegmentRef>()
                    ).add(ref);
                }
            }
        }
    }

    public static final class CombinedResult {
        public final List<ShortReadQuantification> apa;
        public final List<ShortReadIpa.Quantification> ipa;

        CombinedResult(
                List<ShortReadQuantification> apa,
                List<ShortReadIpa.Quantification> ipa
        ) {
            this.apa = apa;
            this.ipa = ipa;
        }
    }

    private static final class IpaWindowRef {
        final int groupIndex;
        final int kind;
        final int start0;
        final int end0;
        final ShortReadIpa.Group group;

        IpaWindowRef(
                int groupIndex,
                int kind,
                ShortReadIpa.Group group,
                CoverageSegment window
        ) {
            this.groupIndex = groupIndex;
            this.kind = kind;
            this.group = group;
            this.start0 = window.start0;
            this.end0 = window.end0;
        }
    }

    private static final class IpaWindowIndex {
        private final Map<String, Map<Integer, List<IpaWindowRef>>> bins =
                new HashMap<String, Map<Integer, List<IpaWindowRef>>>();

        IpaWindowIndex(List<ShortReadIpa.Group> groups) {
            for (int groupIndex = 0; groupIndex < groups.size(); groupIndex++) {
                ShortReadIpa.Group group = groups.get(groupIndex);
                add(groupIndex, group, 0, group.commonWindow);
                add(groupIndex, group, 1, group.ipaSpecificWindow);
                add(groupIndex, group, 2, group.continuationWindow);
            }
        }

        Set<IpaWindowRef> candidates(String chromosome, int start0, int end0) {
            Map<Integer, List<IpaWindowRef>> chromosomeBins = bins.get(chromosome);
            if (chromosomeBins == null) {
                return Collections.emptySet();
            }

            Set<IpaWindowRef> candidates = new HashSet<IpaWindowRef>();
            int firstBin = Math.max(0, start0 / BIN_SIZE);
            int lastBin = Math.max(0, (end0 - 1) / BIN_SIZE);
            for (int bin = firstBin; bin <= lastBin; bin++) {
                List<IpaWindowRef> refs = chromosomeBins.get(bin);
                if (refs != null) {
                    candidates.addAll(refs);
                }
            }
            return candidates;
        }

        private void add(
                int groupIndex,
                ShortReadIpa.Group group,
                int kind,
                CoverageSegment window
        ) {
            IpaWindowRef ref = new IpaWindowRef(groupIndex, kind, group, window);
            int firstBin = window.start0 / BIN_SIZE;
            int lastBin = (window.end0 - 1) / BIN_SIZE;
            Map<Integer, List<IpaWindowRef>> chromosomeBins =
                    bins.computeIfAbsent(
                            group.chrom,
                            key -> new HashMap<Integer, List<IpaWindowRef>>()
                    );
            for (int bin = firstBin; bin <= lastBin; bin++) {
                chromosomeBins.computeIfAbsent(
                        bin,
                        key -> new ArrayList<IpaWindowRef>()
                ).add(ref);
            }
        }
    }

    public static List<ShortReadQuantification> quantify(
            String bam,
            List<ShortReadGroup> groups,
            int minMapq,
            boolean keepDuplicates,
            String library
    ) throws IOException {
        return quantifyWithIpa(
                bam,
                groups,
                Collections.<ShortReadIpa.Group>emptyList(),
                minMapq,
                keepDuplicates,
                library
        ).apa;
    }

    public static CombinedResult quantifyWithIpa(
            String bam,
            List<ShortReadGroup> groups,
            List<ShortReadIpa.Group> ipaGroups,
            int minMapq,
            boolean keepDuplicates,
            String library
    ) throws IOException {
        long[][] alignedBases = initializeAlignedBaseMatrix(groups);
        int[] commonSegmentReads = new int[groups.size()];
        SegmentIndex index = new SegmentIndex(groups);

        int[] ipaCommonReads = new int[ipaGroups.size()];
        int[] ipaSpecificReads = new int[ipaGroups.size()];
        int[] ipaContinuationReads = new int[ipaGroups.size()];
        IpaWindowIndex ipaIndex = new IpaWindowIndex(ipaGroups);

        long alignmentsSeen = 0;
        long alignmentsUsed = 0;
        long alignmentsUsedForIpa = 0;

        try (AlignmentReader reader = AlignmentReaderFactory.open(bam)) {
            AlignmentRecord alignment;
            while ((alignment = reader.next()) != null) {
                alignmentsSeen++;
                if (!isUsable(alignment, minMapq, keepDuplicates)) {
                    continue;
                }

                Strand transcriptStrand = transcriptStrand(alignment, library);
                Set<Integer> commonSegmentHits = new HashSet<Integer>();
                Set<Integer> ipaCommonHits = new HashSet<Integer>();
                Set<Integer> ipaSpecificHits = new HashSet<Integer>();
                Set<Integer> ipaContinuationHits = new HashSet<Integer>();
                boolean overlapsApaRegion = false;

                for (Interval block : alignment.alignedBlocks) {
                    for (SegmentRef segment : index.candidates(
                            alignment.referenceName,
                            block.start0,
                            block.end0
                    )) {
                        if (!strandMatches(
                                segment.group.strand,
                                transcriptStrand,
                                library
                        )) {
                            continue;
                        }

                        int overlapStart = Math.max(block.start0, segment.start0);
                        int overlapEnd = Math.min(block.end0, segment.end0);
                        if (overlapEnd <= overlapStart) {
                            continue;
                        }

                        alignedBases[segment.groupIndex][segment.segmentIndex] +=
                                overlapEnd - overlapStart;
                        if (segment.segmentIndex == 0) {
                            commonSegmentHits.add(segment.groupIndex);
                        }
                        overlapsApaRegion = true;
                    }

                    for (IpaWindowRef window : ipaIndex.candidates(
                            alignment.referenceName,
                            block.start0,
                            block.end0
                    )) {
                        if (!strandMatches(
                                window.group.strand,
                                transcriptStrand,
                                library
                        )) {
                            continue;
                        }

                        int overlapStart = Math.max(block.start0, window.start0);
                        int overlapEnd = Math.min(block.end0, window.end0);
                        if (overlapEnd <= overlapStart) {
                            continue;
                        }

                        if (window.kind == 0) {
                            ipaCommonHits.add(window.groupIndex);
                        } else if (window.kind == 1) {
                            ipaSpecificHits.add(window.groupIndex);
                        } else {
                            ipaContinuationHits.add(window.groupIndex);
                        }
                    }
                }

                for (Integer groupIndex : commonSegmentHits) {
                    commonSegmentReads[groupIndex]++;
                }
                for (Integer groupIndex : ipaCommonHits) {
                    ipaCommonReads[groupIndex]++;
                }
                for (Integer groupIndex : ipaSpecificHits) {
                    ipaSpecificReads[groupIndex]++;
                }
                for (Integer groupIndex : ipaContinuationHits) {
                    ipaContinuationReads[groupIndex]++;
                }

                if (overlapsApaRegion) {
                    alignmentsUsed++;
                }
                if (!ipaCommonHits.isEmpty()
                        || !ipaSpecificHits.isEmpty()
                        || !ipaContinuationHits.isEmpty()) {
                    alignmentsUsedForIpa++;
                }
            }
        }

        System.err.println(
                "[APAURA] short-read alignments seen=" + alignmentsSeen
                        + " overlapping_APA_regions=" + alignmentsUsed
                        + " overlapping_IPA_windows=" + alignmentsUsedForIpa
                        + " file=" + bam
        );

        List<ShortReadQuantification> apaResults =
                new ArrayList<ShortReadQuantification>();
        for (int groupIndex = 0; groupIndex < groups.size(); groupIndex++) {
            apaResults.add(finish(
                    groups.get(groupIndex),
                    alignedBases[groupIndex],
                    commonSegmentReads[groupIndex]
            ));
        }

        List<ShortReadIpa.Quantification> ipaResults =
                new ArrayList<ShortReadIpa.Quantification>();
        for (int groupIndex = 0; groupIndex < ipaGroups.size(); groupIndex++) {
            ipaResults.add(ShortReadIpa.quantify(
                    ipaGroups.get(groupIndex),
                    ipaCommonReads[groupIndex],
                    ipaSpecificReads[groupIndex],
                    ipaContinuationReads[groupIndex]
            ));
        }

        return new CombinedResult(apaResults, ipaResults);
    }

    private static long[][] initializeAlignedBaseMatrix(List<ShortReadGroup> groups) {
        long[][] alignedBases = new long[groups.size()][];
        for (int groupIndex = 0; groupIndex < groups.size(); groupIndex++) {
            alignedBases[groupIndex] =
                    new long[groups.get(groupIndex).segments.size()];
        }
        return alignedBases;
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
        if (!keepDuplicates && alignment.isDuplicate()) {
            return false;
        }
        return !alignment.alignedBlocks.isEmpty();
    }

    private static boolean strandMatches(
            Strand groupStrand,
            Strand transcriptStrand,
            String library
    ) {
        if ("unstranded".equalsIgnoreCase(library)) {
            return true;
        }
        return transcriptStrand == Strand.UNKNOWN || groupStrand == transcriptStrand;
    }

    private static ShortReadQuantification finish(
            ShortReadGroup group,
            long[] alignedBases,
            int commonReads
    ) {
        int segmentCount = group.segments.size();
        double[] rawCoverage = new double[segmentCount];
        double[] weights = new double[segmentCount];

        for (int index = 0; index < segmentCount; index++) {
            int length = group.segments.get(index).length();
            rawCoverage[index] = length > 0
                    ? (double) alignedBases[index] / length
                    : 0.0;
            weights[index] = Math.max(1, length);
        }

        double[] fittedCoverage = IsotonicRegression.decreasing(
                rawCoverage,
                weights
        );
        double[] abundances = inferSiteAbundances(fittedCoverage);
        double[] pau = normalize(abundances);
        int[] pseudoCounts = allocateCounts(pau, commonReads);

        return new ShortReadQuantification(
                group,
                alignedBases,
                rawCoverage,
                fittedCoverage,
                abundances,
                pau,
                pseudoCounts,
                commonReads
        );
    }

    private static double[] inferSiteAbundances(double[] fittedCoverage) {
        double[] abundance = new double[fittedCoverage.length];
        if (fittedCoverage.length == 0) {
            return abundance;
        }

        for (int index = 0; index < fittedCoverage.length - 1; index++) {
            abundance[index] = Math.max(
                    0.0,
                    fittedCoverage[index] - fittedCoverage[index + 1]
            );
        }
        abundance[fittedCoverage.length - 1] = Math.max(
                0.0,
                fittedCoverage[fittedCoverage.length - 1]
        );
        return abundance;
    }

    private static double[] normalize(double[] values) {
        double[] normalized = new double[values.length];
        double total = 0.0;
        for (double value : values) {
            total += value;
        }
        if (total <= 0.0) {
            return normalized;
        }
        for (int index = 0; index < values.length; index++) {
            normalized[index] = values[index] / total;
        }
        return normalized;
    }

    /**
     * Converts fractional PAU estimates to integer effective counts while preserving
     * the total number of reads crossing the common segment.
     */
    private static int[] allocateCounts(double[] probabilities, int totalReads) {
        int[] counts = new int[probabilities.length];
        if (totalReads <= 0 || probabilities.length == 0) {
            return counts;
        }

        double[] fractionalRemainders = new double[probabilities.length];
        int allocated = 0;

        for (int index = 0; index < probabilities.length; index++) {
            double expected = Math.max(0.0, probabilities[index]) * totalReads;
            counts[index] = (int) Math.floor(expected);
            fractionalRemainders[index] = expected - counts[index];
            allocated += counts[index];
        }

        while (allocated < totalReads) {
            int best = indexOfLargest(fractionalRemainders);
            counts[best]++;
            fractionalRemainders[best] = -1.0;
            allocated++;
        }

        while (allocated > totalReads) {
            int best = indexForRemoval(counts, fractionalRemainders);
            if (best < 0) {
                break;
            }
            counts[best]--;
            allocated--;
        }
        return counts;
    }

    private static int indexOfLargest(double[] values) {
        int best = 0;
        for (int index = 1; index < values.length; index++) {
            if (values[index] > values[best]) {
                best = index;
            }
        }
        return best;
    }

    private static int indexForRemoval(int[] counts, double[] remainders) {
        int best = -1;
        for (int index = 0; index < counts.length; index++) {
            if (counts[index] <= 0) {
                continue;
            }
            if (best < 0 || remainders[index] < remainders[best]) {
                best = index;
            }
        }
        return best;
    }

    private static Strand transcriptStrand(AlignmentRecord alignment, String mode) {
        if (mode == null || "unstranded".equalsIgnoreCase(mode)) {
            return Strand.UNKNOWN;
        }

        Strand alignmentStrand = alignment.alignmentStrand();
        if (isForwardMode(mode)) {
            return alignment.isPaired() && alignment.isSecondOfPair()
                    ? alignmentStrand.flip()
                    : alignmentStrand;
        }
        if (isReverseMode(mode)) {
            return alignment.isPaired() && alignment.isSecondOfPair()
                    ? alignmentStrand
                    : alignmentStrand.flip();
        }
        throw new IllegalArgumentException(
                "Unknown --library mode: " + mode
                        + " (use unstranded, forward, or reverse)"
        );
    }

    private static boolean isForwardMode(String mode) {
        return "forward".equalsIgnoreCase(mode)
                || "secondstrand".equalsIgnoreCase(mode)
                || "fr-secondstrand".equalsIgnoreCase(mode);
    }

    private static boolean isReverseMode(String mode) {
        return "reverse".equalsIgnoreCase(mode)
                || "firststrand".equalsIgnoreCase(mode)
                || "fr-firststrand".equalsIgnoreCase(mode);
    }
}
