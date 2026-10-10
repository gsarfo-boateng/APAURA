package org.apaura.analysis;

import org.apaura.model.Sample;
import org.apaura.model.SiteCluster;
import org.apaura.model.Strand;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Applies depth and pooled-usage filters to direct/long-read PAS calls. */
public final class DirectSiteFilter {

    private DirectSiteFilter() {
        // Utility class.
    }

    public static final class Result {
        public final List<SiteCluster> retained;

        /** Gene -> sample -> all annotated end counts before PAS filtering. */
        public final Map<String, Map<String, Integer>> geneTotals;

        Result(
                List<SiteCluster> retained,
                Map<String, Map<String, Integer>> geneTotals
        ) {
            this.retained = retained;
            this.geneTotals = geneTotals;
        }
    }

    public static Result filter(
            List<Sample> samples,
            List<SiteCluster> sites,
            String control,
            String caseCondition,
            int minGeneReads,
            int minSiteReads,
            double minPau
    ) {
        if (minPau < 0 || minPau > 1) {
            throw new IllegalArgumentException(
                    "--min-pau must be a fraction between 0 and 1"
            );
        }

        Map<String, List<SiteCluster>> sitesByGene = groupByGene(sites);
        Map<String, Map<String, Integer>> geneTotals =
                calculateGeneTotals(samples, sitesByGene);

        List<SiteCluster> retained = new ArrayList<SiteCluster>();
        int droppedLowGeneDepth = 0;
        int droppedLowSiteReads = 0;
        int droppedLowPau = 0;

        for (Map.Entry<String, List<SiteCluster>> entry : sitesByGene.entrySet()) {
            Map<String, Integer> totals = geneTotals.get(entry.getKey());
            int controlDepth = pooledConditionDepth(samples, totals, control);
            int caseDepth = pooledConditionDepth(samples, totals, caseCondition);

            if (controlDepth < minGeneReads || caseDepth < minGeneReads) {
                droppedLowGeneDepth += entry.getValue().size();
                continue;
            }

            int pooledGeneDepth = controlDepth + caseDepth;
            for (SiteCluster site : entry.getValue()) {
                int pooledSiteReads = pooledSelectedConditionCount(
                        samples,
                        site,
                        control,
                        caseCondition
                );
                if (pooledSiteReads < minSiteReads) {
                    droppedLowSiteReads++;
                    continue;
                }

                double pooledPau = pooledGeneDepth > 0
                        ? (double) pooledSiteReads / pooledGeneDepth
                        : 0.0;
                if (pooledPau < minPau) {
                    droppedLowPau++;
                    continue;
                }
                retained.add(site);
            }
        }

        reassignRanks(retained);
        System.err.println(
                "[APAURA] PAS filtering retained=" + retained.size()
                        + " dropped_low_gene_depth=" + droppedLowGeneDepth
                        + " dropped_low_site_reads=" + droppedLowSiteReads
                        + " dropped_low_pau=" + droppedLowPau
        );
        return new Result(retained, geneTotals);
    }

    private static Map<String, List<SiteCluster>> groupByGene(List<SiteCluster> sites) {
        Map<String, List<SiteCluster>> grouped =
                new LinkedHashMap<String, List<SiteCluster>>();
        for (SiteCluster site : sites) {
            if (".".equals(site.geneId)) {
                continue;
            }
            grouped.computeIfAbsent(
                    site.geneId,
                    key -> new ArrayList<SiteCluster>()
            ).add(site);
        }
        return grouped;
    }

    private static Map<String, Map<String, Integer>> calculateGeneTotals(
            List<Sample> samples,
            Map<String, List<SiteCluster>> sitesByGene
    ) {
        Map<String, Map<String, Integer>> totals =
                new LinkedHashMap<String, Map<String, Integer>>();

        for (Map.Entry<String, List<SiteCluster>> entry : sitesByGene.entrySet()) {
            Map<String, Integer> bySample = new LinkedHashMap<String, Integer>();
            for (Sample sample : samples) {
                int total = 0;
                for (SiteCluster site : entry.getValue()) {
                    total += count(site, sample.name);
                }
                bySample.put(sample.name, total);
            }
            totals.put(entry.getKey(), bySample);
        }
        return totals;
    }

    private static int pooledConditionDepth(
            List<Sample> samples,
            Map<String, Integer> totals,
            String condition
    ) {
        int depth = 0;
        for (Sample sample : samples) {
            if (sample.condition.equals(condition)) {
                depth += get(totals, sample.name);
            }
        }
        return depth;
    }

    private static int pooledSelectedConditionCount(
            List<Sample> samples,
            SiteCluster site,
            String control,
            String caseCondition
    ) {
        int count = 0;
        for (Sample sample : samples) {
            if (sample.condition.equals(control)
                    || sample.condition.equals(caseCondition)) {
                count += count(site, sample.name);
            }
        }
        return count;
    }

    private static int count(SiteCluster site, String sample) {
        Integer value = site.counts.get(sample);
        return value == null ? 0 : value;
    }

    private static int get(Map<String, Integer> values, String key) {
        Integer value = values.get(key);
        return value == null ? 0 : value;
    }

    private static void reassignRanks(List<SiteCluster> sites) {
        Map<String, List<SiteCluster>> byGene =
                new LinkedHashMap<String, List<SiteCluster>>();
        for (SiteCluster site : sites) {
            byGene.computeIfAbsent(
                    site.geneId,
                    key -> new ArrayList<SiteCluster>()
            ).add(site);
        }

        for (List<SiteCluster> geneSites : byGene.values()) {
            geneSites.sort((left, right) -> {
                int comparison = Integer.compare(
                        left.representativePos0,
                        right.representativePos0
                );
                return left.strand == Strand.MINUS ? -comparison : comparison;
            });

            int siteCount = geneSites.size();
            for (int index = 0; index < siteCount; index++) {
                SiteCluster site = geneSites.get(index);
                site.siteRank = index + 1;
                site.normalizedPosition = siteCount == 1
                        ? 0.5
                        : (double) index / (siteCount - 1);
            }
        }
    }
}
