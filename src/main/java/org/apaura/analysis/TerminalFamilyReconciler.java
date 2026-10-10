package org.apaura.analysis;

import org.apaura.model.Strand;
import org.apaura.model.TranscriptModel;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Reconciles transcript annotations into non-redundant terminal-exon families.
 *
 * <p>Annotated transcript ends are first clustered gene-wide into physical PAS
 * clusters. Terminal-anchor buckets that share a physical PAS are then merged
 * into one family. This prevents the same genomic PAS from being quantified as
 * multiple independent tandem-UTR events simply because transcript annotations
 * use different terminal-exon 5' boundaries.</p>
 */
final class TerminalFamilyReconciler {

    private TerminalFamilyReconciler() {
        // Utility class.
    }

    static final class Family {
        final String geneId;
        final String geneName;
        final String biotype;
        final String chrom;
        final Strand strand;
        final int anchor0;
        final boolean coding;
        final int sourceAnchorBuckets;
        final int transcriptCount;
        final List<Integer> sourceAnchors;
        final List<Integer> pasPositions;

        Family(
                String geneId,
                String geneName,
                String biotype,
                String chrom,
                Strand strand,
                int anchor0,
                boolean coding,
                int sourceAnchorBuckets,
                int transcriptCount,
                List<Integer> sourceAnchors,
                List<Integer> pasPositions
        ) {
            this.geneId = geneId;
            this.geneName = geneName;
            this.biotype = biotype;
            this.chrom = chrom;
            this.strand = strand;
            this.anchor0 = anchor0;
            this.coding = coding;
            this.sourceAnchorBuckets = sourceAnchorBuckets;
            this.transcriptCount = transcriptCount;
            this.sourceAnchors = Collections.unmodifiableList(
                    new ArrayList<Integer>(sourceAnchors)
            );
            this.pasPositions = Collections.unmodifiableList(
                    new ArrayList<Integer>(pasPositions)
            );
        }
    }

    static final class Result {
        final List<Family> families;
        final int initialAnchorBuckets;
        final int mergedAnchorBuckets;

        Result(
                List<Family> families,
                int initialAnchorBuckets,
                int mergedAnchorBuckets
        ) {
            this.families = families;
            this.initialAnchorBuckets = initialAnchorBuckets;
            this.mergedAnchorBuckets = mergedAnchorBuckets;
        }
    }

    private static final class AnchorBucket {
        String geneId;
        String geneName;
        String biotype;
        String chrom;
        Strand strand;
        int representativeAnchor0;
        int conservativeAnchor0;
        boolean coding;
        final List<TranscriptModel> transcripts = new ArrayList<TranscriptModel>();
        final Set<Integer> canonicalPas = new LinkedHashSet<Integer>();

        AnchorBucket(TranscriptModel transcript) {
            geneId = transcript.geneId;
            geneName = transcript.geneName;
            biotype = transcript.biotype;
            chrom = transcript.chrom;
            strand = transcript.strand;
            representativeAnchor0 = transcript.terminalExonAnchor0();
            conservativeAnchor0 = transcript.terminalExonAnchor0();
            coding = transcript.hasCds();
            transcripts.add(transcript);
        }
    }

    static Result reconcile(
            List<TranscriptModel> geneTranscripts,
            int anchorTolerance,
            int pasMergeDistance
    ) {
        if (geneTranscripts.isEmpty()) {
            return new Result(Collections.<Family>emptyList(), 0, 0);
        }

        List<TranscriptModel> ordered = new ArrayList<TranscriptModel>(geneTranscripts);
        ordered.sort(Comparator.comparingInt(TranscriptModel::terminalExonAnchor0));

        List<AnchorBucket> buckets = bucketByTerminalAnchor(
                ordered,
                anchorTolerance
        );
        Map<String, Map<Integer, Integer>> canonicalByLocus =
                canonicalPasByLocus(ordered, pasMergeDistance);
        assignCanonicalPas(buckets, canonicalByLocus);

        DisjointSet disjointSet = new DisjointSet(buckets.size());
        Map<String, Integer> firstBucketByPas = new HashMap<String, Integer>();
        for (int bucketIndex = 0; bucketIndex < buckets.size(); bucketIndex++) {
            AnchorBucket bucket = buckets.get(bucketIndex);
            for (Integer pas : bucket.canonicalPas) {
                String key = locusKey(bucket.chrom, bucket.strand) + "\u0001" + pas;
                Integer previous = firstBucketByPas.putIfAbsent(key, bucketIndex);
                if (previous != null) {
                    disjointSet.union(previous, bucketIndex);
                }
            }
        }

        Map<Integer, List<AnchorBucket>> components =
                new LinkedHashMap<Integer, List<AnchorBucket>>();
        for (int bucketIndex = 0; bucketIndex < buckets.size(); bucketIndex++) {
            int root = disjointSet.find(bucketIndex);
            components.computeIfAbsent(
                    root,
                    key -> new ArrayList<AnchorBucket>()
            ).add(buckets.get(bucketIndex));
        }

        List<Family> families = new ArrayList<Family>();
        for (List<AnchorBucket> component : components.values()) {
            families.add(toFamily(component));
        }

        int merged = buckets.size() - families.size();
        return new Result(families, buckets.size(), merged);
    }

    private static List<AnchorBucket> bucketByTerminalAnchor(
            List<TranscriptModel> transcripts,
            int tolerance
    ) {
        List<AnchorBucket> buckets = new ArrayList<AnchorBucket>();
        for (TranscriptModel transcript : transcripts) {
            AnchorBucket bucket = nearestBucket(transcript, buckets, tolerance);
            if (bucket == null) {
                buckets.add(new AnchorBucket(transcript));
            } else {
                mergeTranscriptIntoBucket(bucket, transcript);
            }
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
                    bucket.representativeAnchor0 - transcript.terminalExonAnchor0()
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
        int transcriptAnchor = transcript.terminalExonAnchor0();
        bucket.representativeAnchor0 =
                (bucket.representativeAnchor0 + transcriptAnchor) / 2;
        if (bucket.strand == Strand.MINUS) {
            bucket.conservativeAnchor0 = Math.min(
                    bucket.conservativeAnchor0,
                    transcriptAnchor
            );
        } else {
            bucket.conservativeAnchor0 = Math.max(
                    bucket.conservativeAnchor0,
                    transcriptAnchor
            );
        }
        bucket.coding |= transcript.hasCds();
        bucket.biotype = preferredBiotype(bucket.biotype, transcript.biotype);
        bucket.transcripts.add(transcript);
    }

    private static Map<String, Map<Integer, Integer>> canonicalPasByLocus(
            List<TranscriptModel> transcripts,
            int mergeDistance
    ) {
        Map<String, List<Integer>> positionsByLocus =
                new LinkedHashMap<String, List<Integer>>();
        for (TranscriptModel transcript : transcripts) {
            positionsByLocus.computeIfAbsent(
                    locusKey(transcript.chrom, transcript.strand),
                    key -> new ArrayList<Integer>()
            ).add(transcript.pasPos0());
        }

        Map<String, Map<Integer, Integer>> result =
                new HashMap<String, Map<Integer, Integer>>();
        for (Map.Entry<String, List<Integer>> entry : positionsByLocus.entrySet()) {
            result.put(entry.getKey(), canonicalize(entry.getValue(), mergeDistance));
        }
        return result;
    }

    private static Map<Integer, Integer> canonicalize(
            List<Integer> positions,
            int mergeDistance
    ) {
        List<Integer> unique = new ArrayList<Integer>(
                new LinkedHashSet<Integer>(positions)
        );
        Collections.sort(unique);

        Map<Integer, Integer> canonical = new HashMap<Integer, Integer>();
        if (unique.isEmpty()) {
            return canonical;
        }

        List<Integer> cluster = new ArrayList<Integer>();
        cluster.add(unique.get(0));
        for (int index = 1; index < unique.size(); index++) {
            int position = unique.get(index);
            if (position - unique.get(index - 1) <= mergeDistance) {
                cluster.add(position);
            } else {
                assignCluster(cluster, canonical);
                cluster.clear();
                cluster.add(position);
            }
        }
        assignCluster(cluster, canonical);
        return canonical;
    }

    private static void assignCluster(
            List<Integer> cluster,
            Map<Integer, Integer> canonical
    ) {
        int representative = cluster.get(cluster.size() / 2);
        for (Integer position : cluster) {
            canonical.put(position, representative);
        }
    }

    private static void assignCanonicalPas(
            List<AnchorBucket> buckets,
            Map<String, Map<Integer, Integer>> canonicalByLocus
    ) {
        for (AnchorBucket bucket : buckets) {
            Map<Integer, Integer> canonical = canonicalByLocus.get(
                    locusKey(bucket.chrom, bucket.strand)
            );
            for (TranscriptModel transcript : bucket.transcripts) {
                Integer position = canonical.get(transcript.pasPos0());
                if (position != null) {
                    bucket.canonicalPas.add(position);
                }
            }
        }
    }

        private static Family toFamily(List<AnchorBucket> buckets) {
            AnchorBucket first = buckets.get(0);

            Set<Integer> positions = new LinkedHashSet<Integer>();
            Set<Integer> sourceAnchors = new LinkedHashSet<Integer>();

            boolean coding = false;
            String biotype = ".";
            int transcriptCount = 0;
            int anchor0 = first.conservativeAnchor0;

            for (AnchorBucket bucket : buckets) {
                positions.addAll(bucket.canonicalPas);
                coding |= bucket.coding;
                biotype = preferredBiotype(biotype, bucket.biotype);
                transcriptCount += bucket.transcripts.size();

                for (TranscriptModel transcript : bucket.transcripts) {
                    sourceAnchors.add(transcript.terminalExonAnchor0());
                }

                if (buckets.size() > 1) {
                    if (first.strand == Strand.MINUS) {
                        anchor0 = Math.min(
                                anchor0,
                                bucket.conservativeAnchor0
                        );
                    } else {
                        anchor0 = Math.max(
                                anchor0,
                                bucket.conservativeAnchor0
                        );
                    }
                }
            }

            List<Integer> orderedSourceAnchors =
                    new ArrayList<Integer>(sourceAnchors);
            Collections.sort(orderedSourceAnchors);

            List<Integer> orderedPositions =
                    new ArrayList<Integer>(positions);
            Collections.sort(orderedPositions);

            if (first.strand == Strand.MINUS) {
                Collections.reverse(orderedPositions);
            }

            return new Family(
                    first.geneId,
                    first.geneName,
                    biotype,
                    first.chrom,
                    first.strand,
                    anchor0,
                    coding,
                    buckets.size(),
                    transcriptCount,
                    orderedSourceAnchors,
                    orderedPositions
            );
        }

    private static String preferredBiotype(String left, String right) {
        if ("protein_coding".equals(left) || "protein_coding".equals(right)) {
            return "protein_coding";
        }
        if (left == null || ".".equals(left)) {
            return right == null ? "." : right;
        }
        return left;
    }

    private static String locusKey(String chrom, Strand strand) {
        return chrom + "\u0001" + strand.symbol();
    }

    private static final class DisjointSet {
        private final int[] parent;
        private final int[] rank;

        DisjointSet(int size) {
            parent = new int[size];
            rank = new int[size];
            for (int index = 0; index < size; index++) {
                parent[index] = index;
            }
        }

        int find(int value) {
            if (parent[value] != value) {
                parent[value] = find(parent[value]);
            }
            return parent[value];
        }

        void union(int left, int right) {
            int leftRoot = find(left);
            int rightRoot = find(right);
            if (leftRoot == rightRoot) {
                return;
            }

            if (rank[leftRoot] < rank[rightRoot]) {
                parent[leftRoot] = rightRoot;
            } else if (rank[leftRoot] > rank[rightRoot]) {
                parent[rightRoot] = leftRoot;
            } else {
                parent[rightRoot] = leftRoot;
                rank[leftRoot]++;
            }
        }
    }
}