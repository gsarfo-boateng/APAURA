package org.apaura.analysis;

import org.apaura.model.EndKey;
import org.apaura.model.Sample;
import org.apaura.model.SiteCluster;
import org.apaura.model.Strand;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/** Clusters nearby transcript-end coordinates into candidate poly(A) sites. */
public final class SiteClusterer {
    private SiteClusterer() {
    }

    private static final class PoolPoint {
        final EndKey key;
        int count;

        PoolPoint(EndKey key, int count) {
            this.key = key;
            this.count = count;
        }
    }

    public static List<SiteCluster> cluster(
            List<Sample> samples,
            Map<String, Map<EndKey, Integer>> perSample,
            int distance
    ) {
        return cluster(samples, perSample, distance, "mode");
    }

    /**
     * Clusters coordinates separately by chromosome and strand.
     * Representative modes are "mode" (highest pooled support) and
     * "median" (weighted median coordinate).
     */
    public static List<SiteCluster> cluster(
            List<Sample> samples,
            Map<String, Map<EndKey, Integer>> perSample,
            int distance,
            String representativeMode
    ) {
        List<PoolPoint> pooledPoints = poolCoordinates(samples, perSample);
        Collections.sort(pooledPoints, new Comparator<PoolPoint>() {
            @Override
            public int compare(PoolPoint left, PoolPoint right) {
                return left.key.compareTo(right.key);
            }
        });

        List<SiteCluster> clusters = new ArrayList<SiteCluster>();
        int serial = 1;
        int startIndex = 0;

        while (startIndex < pooledPoints.size()) {
            int endIndex = findClusterEnd(pooledPoints, startIndex, distance);
            List<PoolPoint> group = pooledPoints.subList(startIndex, endIndex);

            int representative = chooseRepresentative(group, representativeMode);
            PoolPoint first = group.get(0);
            PoolPoint last = group.get(group.size() - 1);

            SiteCluster cluster = new SiteCluster(
                    String.format(Locale.ROOT, "PAS%07d", serial++),
                    first.key.chrom,
                    first.key.strand,
                    first.key.pos0,
                    last.key.pos0,
                    representative
            );

            addSampleCounts(cluster, samples, perSample, group);
            clusters.add(cluster);
            startIndex = endIndex;
        }

        return clusters;
    }

    private static List<PoolPoint> poolCoordinates(
            List<Sample> samples,
            Map<String, Map<EndKey, Integer>> perSample
    ) {
        Map<EndKey, Integer> pooledCounts = new HashMap<EndKey, Integer>();

        for (Sample sample : samples) {
            Map<EndKey, Integer> sampleCounts = perSample.get(sample.name);
            if (sampleCounts == null) {
                continue;
            }

            for (Map.Entry<EndKey, Integer> entry : sampleCounts.entrySet()) {
                Integer previous = pooledCounts.get(entry.getKey());
                pooledCounts.put(
                        entry.getKey(),
                        (previous == null ? 0 : previous) + entry.getValue()
                );
            }
        }

        List<PoolPoint> points = new ArrayList<PoolPoint>();
        for (Map.Entry<EndKey, Integer> entry : pooledCounts.entrySet()) {
            points.add(new PoolPoint(entry.getKey(), entry.getValue()));
        }
        return points;
    }

    private static int findClusterEnd(
            List<PoolPoint> points,
            int startIndex,
            int distance
    ) {
        PoolPoint first = points.get(startIndex);
        String chromosome = first.key.chrom;
        Strand strand = first.key.strand;
        int previousPosition = first.key.pos0;
        int endIndex = startIndex + 1;

        while (endIndex < points.size()) {
            PoolPoint current = points.get(endIndex);
            if (!current.key.chrom.equals(chromosome)
                    || current.key.strand != strand
                    || current.key.pos0 - previousPosition > distance) {
                break;
            }

            previousPosition = current.key.pos0;
            endIndex++;
        }
        return endIndex;
    }

    private static int chooseRepresentative(
            List<PoolPoint> group,
            String representativeMode
    ) {
        if ("median".equalsIgnoreCase(representativeMode)
                || "weighted-median".equalsIgnoreCase(representativeMode)) {
            return weightedMedian(group);
        }
        if ("mode".equalsIgnoreCase(representativeMode)
                || "peak".equalsIgnoreCase(representativeMode)) {
            return mode(group);
        }
        throw new IllegalArgumentException(
                "Unknown cluster representative mode: "
                        + representativeMode
                        + " (use mode or median)"
        );
    }

    private static void addSampleCounts(
            SiteCluster cluster,
            List<Sample> samples,
            Map<String, Map<EndKey, Integer>> perSample,
            List<PoolPoint> group
    ) {
        for (Sample sample : samples) {
            int count = 0;
            Map<EndKey, Integer> sampleCounts = perSample.get(sample.name);
            if (sampleCounts != null) {
                for (PoolPoint point : group) {
                    Integer value = sampleCounts.get(point.key);
                    if (value != null) {
                        count += value;
                    }
                }
            }
            cluster.counts.put(sample.name, count);
        }
    }

    private static int weightedMedian(List<PoolPoint> group) {
        int total = 0;
        for (PoolPoint point : group) {
            total += point.count;
        }

        int halfway = (total + 1) / 2;
        int cumulative = 0;
        for (PoolPoint point : group) {
            cumulative += point.count;
            if (cumulative >= halfway) {
                return point.key.pos0;
            }
        }
        return group.get(0).key.pos0;
    }

    private static int mode(List<PoolPoint> group) {
        PoolPoint best = group.get(0);
        for (PoolPoint point : group) {
            if (point.count > best.count
                    || (point.count == best.count && point.key.pos0 < best.key.pos0)) {
                best = point;
            }
        }
        return best.key.pos0;
    }
}
