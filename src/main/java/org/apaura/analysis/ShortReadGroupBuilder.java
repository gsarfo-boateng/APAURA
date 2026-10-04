package org.apaura.analysis;

import org.apaura.model.CoverageSegment;
import org.apaura.model.ShortReadGroup;
import org.apaura.model.ShortReadSite;
import org.apaura.model.Strand;
import org.apaura.model.TranscriptModel;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/** Builds non-redundant tandem terminal-exon APA groups from transcript annotations. */
public final class ShortReadGroupBuilder {

    private ShortReadGroupBuilder() {
        // Utility class.
    }

    public static final class BuildResult {
        public final List<ShortReadGroup> groups;
        public final int initialAnchorGroups;
        public final int reconciledFamilies;
        public final int mergedAnchorGroups;
        public final int skippedNonQuantifiableFamilies;

        BuildResult(
                List<ShortReadGroup> groups,
                int initialAnchorGroups,
                int reconciledFamilies,
                int mergedAnchorGroups,
                int skippedNonQuantifiableFamilies
        ) {
            this.groups = groups;
            this.initialAnchorGroups = initialAnchorGroups;
            this.reconciledFamilies = reconciledFamilies;
            this.mergedAnchorGroups = mergedAnchorGroups;
            this.skippedNonQuantifiableFamilies = skippedNonQuantifiableFamilies;
        }
    }

    public static List<ShortReadGroup> build(
            List<TranscriptModel> transcripts,
            int anchorTolerance,
            int pasMergeDistance,
            int minSegmentLength
    ) {
        return buildWithReport(
                transcripts,
                anchorTolerance,
                pasMergeDistance,
                minSegmentLength
        ).groups;
    }

    public static BuildResult buildWithReport(
            List<TranscriptModel> transcripts,
            int anchorTolerance,
            int pasMergeDistance,
            int minSegmentLength
    ) {
        Map<String, List<TranscriptModel>> transcriptsByGene = groupByGene(transcripts);
        List<ShortReadGroup> groups = new ArrayList<ShortReadGroup>();
        int initialAnchorGroups = 0;
        int reconciledFamilies = 0;
        int mergedAnchorGroups = 0;
        int skippedNonQuantifiableFamilies = 0;

        for (List<TranscriptModel> geneTranscripts : transcriptsByGene.values()) {
            TerminalFamilyReconciler.Result reconciliation =
                    TerminalFamilyReconciler.reconcile(
                            geneTranscripts,
                            anchorTolerance,
                            pasMergeDistance
                    );
            initialAnchorGroups += reconciliation.initialAnchorBuckets;
            reconciledFamilies += reconciliation.families.size();
            mergedAnchorGroups += reconciliation.mergedAnchorBuckets;
            skippedNonQuantifiableFamilies += addGroups(
                    reconciliation.families,
                    minSegmentLength,
                    groups
            );
        }

        groups.sort((left, right) -> {
            int comparison = left.chrom.compareTo(right.chrom);
            if (comparison != 0) {
                return comparison;
            }
            comparison = Integer.compare(left.anchor0, right.anchor0);
            return comparison != 0
                    ? comparison
                    : left.id.compareTo(right.id);
        });
        assertUniquePasOwnership(groups);

        return new BuildResult(
                Collections.unmodifiableList(groups),
                initialAnchorGroups,
                reconciledFamilies,
                mergedAnchorGroups,
                skippedNonQuantifiableFamilies
        );
    }

    private static Map<String, List<TranscriptModel>> groupByGene(
            List<TranscriptModel> transcripts
    ) {
        Map<String, List<TranscriptModel>> grouped =
                new LinkedHashMap<String, List<TranscriptModel>>();
        for (TranscriptModel transcript : transcripts) {
            grouped.computeIfAbsent(
                    transcript.geneId,
                    key -> new ArrayList<TranscriptModel>()
            ).add(transcript);
        }
        return grouped;
    }

    private static int addGroups(
            List<TerminalFamilyReconciler.Family> families,
            int minSegmentLength,
            List<ShortReadGroup> output
    ) {
        int groupNumber = 0;
        int skipped = 0;

        for (TerminalFamilyReconciler.Family family : families) {
            if (family.pasPositions.size() < 2) {
                continue;
            }

            List<CoverageSegment> coverageSegments = segments(
                    family.anchor0,
                    family.pasPositions,
                    family.strand
            );
            if (!validSegments(
                    coverageSegments,
                    family.pasPositions.size(),
                    minSegmentLength
            )) {
                skipped++;
                continue;
            }

            groupNumber++;
            output.add(createGroup(
                    family,
                    coverageSegments,
                    groupNumber
            ));
        }
        return skipped;
    }

    private static boolean validSegments(
            List<CoverageSegment> segments,
            int expectedCount,
            int minSegmentLength
    ) {
        if (segments.size() != expectedCount) {
            return false;
        }
        for (CoverageSegment segment : segments) {
            if (segment.length() < minSegmentLength) {
                return false;
            }
        }
        return true;
    }

    private static ShortReadGroup createGroup(
            TerminalFamilyReconciler.Family family,
            List<CoverageSegment> segments,
            int groupNumber
    ) {
        String id = family.geneId
                + "|TE"
                + String.format(Locale.ROOT, "%03d", groupNumber);
        List<ShortReadSite> sites = new ArrayList<ShortReadSite>();
        int siteCount = family.pasPositions.size();

        for (int index = 0; index < siteCount; index++) {
            sites.add(new ShortReadSite(
                    id + "|PAS" + (index + 1),
                    family.pasPositions.get(index),
                    index + 1,
                    siteCount == 1 ? 0.5 : (double) index / (siteCount - 1)
            ));
        }

        String scope = family.coding
                ? "TANDEM_3UTR"
                : "TANDEM_TERMINAL_EXON_NONCODING";
        return new ShortReadGroup(
                id,
                family.geneId,
                family.geneName,
                family.biotype,
                family.chrom,
                family.strand,
                scope,
                family.anchor0,
                sites,
                segments
        );
    }

    private static void assertUniquePasOwnership(List<ShortReadGroup> groups) {
        Map<String, String> ownerByPas = new LinkedHashMap<String, String>();
        for (ShortReadGroup group : groups) {
            for (ShortReadSite site : group.sites) {
                String key = group.geneId
                        + "\u0001"
                        + group.chrom
                        + "\u0001"
                        + group.strand.symbol()
                        + "\u0001"
                        + site.pos0;
                String previous = ownerByPas.putIfAbsent(key, group.id);
                if (previous != null && !previous.equals(group.id)) {
                    throw new IllegalStateException(
                            "Internal terminal-family error: PAS "
                                    + group.chrom
                                    + ":"
                                    + (site.pos0 + 1)
                                    + group.strand.symbol()
                                    + " belongs to both "
                                    + previous
                                    + " and "
                                    + group.id
                    );
                }
            }
        }
    }

    private static List<CoverageSegment> segments(
            int anchor0,
            List<Integer> pasPositions,
            Strand strand
    ) {
        return strand == Strand.MINUS
                ? minusStrandSegments(anchor0, pasPositions)
                : plusStrandSegments(anchor0, pasPositions);
    }

    private static List<CoverageSegment> minusStrandSegments(
            int anchor0,
            List<Integer> pasPositions
    ) {
        List<CoverageSegment> segments = new ArrayList<CoverageSegment>();
        int firstPas = pasPositions.get(0);
        if (anchor0 <= firstPas) {
            return segments;
        }

        segments.add(new CoverageSegment(0, firstPas, anchor0, "COMMON"));
        for (int index = 1; index < pasPositions.size(); index++) {
            int start0 = pasPositions.get(index);
            int end0 = pasPositions.get(index - 1);
            if (end0 <= start0) {
                return Collections.emptyList();
            }
            segments.add(new CoverageSegment(
                    index,
                    start0,
                    end0,
                    "DISTAL_" + index
            ));
        }
        return segments;
    }

    private static List<CoverageSegment> plusStrandSegments(
            int anchor0,
            List<Integer> pasPositions
    ) {
        List<CoverageSegment> segments = new ArrayList<CoverageSegment>();
        int firstPas = pasPositions.get(0);
        if (firstPas + 1 <= anchor0) {
            return segments;
        }

        segments.add(new CoverageSegment(0, anchor0, firstPas + 1, "COMMON"));
        for (int index = 1; index < pasPositions.size(); index++) {
            int start0 = pasPositions.get(index - 1) + 1;
            int end0 = pasPositions.get(index) + 1;
            if (end0 <= start0) {
                return Collections.emptyList();
            }
            segments.add(new CoverageSegment(
                    index,
                    start0,
                    end0,
                    "DISTAL_" + index
            ));
        }
        return segments;
    }
}