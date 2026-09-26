package org.apaura.analysis;

import org.apaura.model.CoverageSegment;
import org.apaura.model.ShortReadGroup;
import org.apaura.model.ShortReadSite;
import org.apaura.model.Strand;
import org.apaura.model.TranscriptModel;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/** Builds tandem terminal-exon APA groups from transcript annotations. */
public final class ShortReadGroupBuilder {

    private ShortReadGroupBuilder() {
        // Utility class.
    }

    private static final class AnchorBucket {
        String geneId;
        String geneName;
        String biotype;
        String chrom;
        Strand strand;
        int anchor0;
        boolean coding;
        final List<Integer> pasPositions = new ArrayList<Integer>();

        AnchorBucket(TranscriptModel transcript) {
            geneId = transcript.geneId;
            geneName = transcript.geneName;
            biotype = transcript.biotype;
            chrom = transcript.chrom;
            strand = transcript.strand;
            anchor0 = transcript.terminalExonAnchor0();
            coding = transcript.hasCds();
        }
    }

    public static List<ShortReadGroup> build(
            List<TranscriptModel> transcripts,
            int anchorTolerance,
            int pasMergeDistance,
            int minSegmentLength
    ) {
        Map<String, List<TranscriptModel>> transcriptsByGene = groupByGene(transcripts);
        List<ShortReadGroup> groups = new ArrayList<ShortReadGroup>();

        for (List<TranscriptModel> geneTranscripts : transcriptsByGene.values()) {
            geneTranscripts.sort(Comparator.comparingInt(
                    TranscriptModel::terminalExonAnchor0
            ));
            List<AnchorBucket> buckets = bucketByTerminalAnchor(
                    geneTranscripts,
                    anchorTolerance
            );
            addGroups(
                    buckets,
                    pasMergeDistance,
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
        return groups;
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

    private static List<AnchorBucket> bucketByTerminalAnchor(
            List<TranscriptModel> transcripts,
            int tolerance
    ) {
        List<AnchorBucket> buckets = new ArrayList<AnchorBucket>();
        for (TranscriptModel transcript : transcripts) {
            AnchorBucket bucket = nearestBucket(transcript, buckets, tolerance);
            if (bucket == null) {
                bucket = new AnchorBucket(transcript);
                buckets.add(bucket);
            } else {
                mergeTranscriptIntoBucket(bucket, transcript);
            }
            bucket.pasPositions.add(transcript.pasPos0());
        }
        return buckets;
    }

    private static AnchorBucket nearestBucket(
            TranscriptModel transcript,
            List<AnchorBucket> buckets,
            int tolerance
    ) {
        AnchorBucket best = null;
        int bestDistance = Integer.MAX_VALUE;

        for (AnchorBucket bucket : buckets) {
            if (!bucket.chrom.equals(transcript.chrom)
                    || bucket.strand != transcript.strand) {
                continue;
            }
            int distance = Math.abs(
                    bucket.anchor0 - transcript.terminalExonAnchor0()
            );
            if (distance <= tolerance && distance < bestDistance) {
                best = bucket;
                bestDistance = distance;
            }
        }
        return best;
    }

    private static void mergeTranscriptIntoBucket(
            AnchorBucket bucket,
            TranscriptModel transcript
    ) {
        bucket.anchor0 = (bucket.anchor0 + transcript.terminalExonAnchor0()) / 2;
        bucket.coding |= transcript.hasCds();
        if (".".equals(bucket.biotype) && !".".equals(transcript.biotype)) {
            bucket.biotype = transcript.biotype;
        }
    }

    private static void addGroups(
            List<AnchorBucket> buckets,
            int pasMergeDistance,
            int minSegmentLength,
            List<ShortReadGroup> output
    ) {
        int groupNumber = 0;
        for (AnchorBucket bucket : buckets) {
            List<Integer> pasPositions = mergePositions(
                    bucket.pasPositions,
                    bucket.strand,
                    pasMergeDistance
            );
            if (pasPositions.size() < 2) {
                continue;
            }

            List<CoverageSegment> segments = segments(
                    bucket.anchor0,
                    pasPositions,
                    bucket.strand
            );
            if (!validSegments(segments, pasPositions.size(), minSegmentLength)) {
                continue;
            }

            groupNumber++;
            output.add(createGroup(bucket, pasPositions, segments, groupNumber));
        }
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
            AnchorBucket bucket,
            List<Integer> pasPositions,
            List<CoverageSegment> segments,
            int groupNumber
    ) {
        String id = bucket.geneId
                + "|TE"
                + String.format(Locale.ROOT, "%03d", groupNumber);
        List<ShortReadSite> sites = new ArrayList<ShortReadSite>();
        int siteCount = pasPositions.size();

        for (int index = 0; index < siteCount; index++) {
            sites.add(new ShortReadSite(
                    id + "|PAS" + (index + 1),
                    pasPositions.get(index),
                    index + 1,
                    siteCount == 1 ? 0.5 : (double) index / (siteCount - 1)
            ));
        }

        String scope = bucket.coding
                ? "TANDEM_3UTR"
                : "TANDEM_TERMINAL_EXON_NONCODING";
        return new ShortReadGroup(
                id,
                bucket.geneId,
                bucket.geneName,
                bucket.biotype,
                bucket.chrom,
                bucket.strand,
                scope,
                bucket.anchor0,
                sites,
                segments
        );
    }

    private static List<Integer> mergePositions(
            List<Integer> positions,
            Strand strand,
            int distance
    ) {
        if (positions.isEmpty()) {
            return Collections.emptyList();
        }

        List<Integer> unique = new ArrayList<Integer>(
                new LinkedHashSet<Integer>(positions)
        );
        Collections.sort(unique);

        List<List<Integer>> clusters = new ArrayList<List<Integer>>();
        List<Integer> current = new ArrayList<Integer>();
        current.add(unique.get(0));

        for (int index = 1; index < unique.size(); index++) {
            if (unique.get(index) - unique.get(index - 1) <= distance) {
                current.add(unique.get(index));
            } else {
                clusters.add(current);
                current = new ArrayList<Integer>();
                current.add(unique.get(index));
            }
        }
        clusters.add(current);

        List<Integer> merged = new ArrayList<Integer>();
        for (List<Integer> cluster : clusters) {
            Collections.sort(cluster);
            merged.add(cluster.get(cluster.size() / 2));
        }
        Collections.sort(merged);
        if (strand == Strand.MINUS) {
            Collections.reverse(merged);
        }
        return merged;
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
