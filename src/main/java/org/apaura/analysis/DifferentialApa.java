package org.apaura.analysis;

import org.apaura.model.Sample;
import org.apaura.model.SiteCluster;
import org.apaura.model.Strand;
import org.apaura.stats.BetaBinomial;
import org.apaura.stats.FisherExact;
import org.apaura.stats.MultipleTesting;
import org.apaura.stats.TieAdjustedKs;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/** Performs site-level and gene-level differential APA analysis. */
public final class DifferentialApa {

    private DifferentialApa() {
        // Utility class.
    }

    public static final class SiteResult {
        public SiteCluster site;
        public double controlPau;
        public double casePau;
        public double delta;
        public double p;
        public double q;
        public double rho;
        public String method;
    }

    public static final class GeneResult {
        public String geneId;
        public String geneName;
        public String biotype;
        public String apaType;
        public String apaDirection = "stable";
        public String apaEvent = "none:stable";
        public boolean apaChange = false;
        public String tutrDirection = "not_applicable";
        public String aleDirection = "not_applicable";
        public String ipaDirection = "not_applicable";
        public String internalExonDirection = "not_applicable";
        public String architectureComponents = "none";
        public String changedMechanisms = "none";
        public String mechanismDirections = "none";
        public String dominantControlArchitecture = "not_applicable";
        public String dominantCaseArchitecture = "not_applicable";
        public String architectureTransition = "not_applicable";
        public int nSites;
        public int ksMaxDiffAt0;
        public int nTerminalBranches;
        public int nTutrBranches;
        public double controlWpi;
        public double caseWpi;
        public double deltaWpi;
        public double apd;
        public double minFdr;
        public double ksD;
        public double ksSignedD;
        public double ksZ;
        public double ksPvalue;
        public double ksFdr;
        public double controlTerminalUsage;
        public double caseTerminalUsage;
        public double deltaTerminalUsage;
        public double controlIpaUsage;
        public double caseIpaUsage;
        public double deltaIpaUsage;
        public double terminationShiftIndex;
        public double terminationJsd;
        public double aleSwitchIndex;
        public double aleJsd;
        public double controlAleWpi = Double.NaN;
        public double caseAleWpi = Double.NaN;
        public double deltaAleWpi = Double.NaN;
        public double controlTutrWpi = Double.NaN;
        public double caseTutrWpi = Double.NaN;
        public double deltaTutrWpi = Double.NaN;
        public double tutrShiftMagnitude;
    }

    public static final class Results {
        public List<SiteResult> sites = new ArrayList<SiteResult>();
        public List<GeneResult> genes = new ArrayList<GeneResult>();
    }

    public static Results run(
            List<Sample> samples,
            List<SiteCluster> sites,
            String control,
            String caseCondition,
            int minTotal,
            double effect,
            double fdr,
            double wpiThreshold
    ) {
        return run(
                samples,
                sites,
                control,
                caseCondition,
                minTotal,
                effect,
                fdr,
                wpiThreshold,
                null
        );
    }

    /**
     * Runs differential APA analysis.
     *
     * @param geneTotals optional gene -> sample -> total annotated end counts before PAS filtering
     */
    public static Results run(
            List<Sample> samples,
            List<SiteCluster> sites,
            String control,
            String caseCondition,
            int minTotal,
            double effect,
            double fdr,
            double wpiThreshold,
            Map<String, Map<String, Integer>> geneTotals
    ) {
        Map<String, List<SiteCluster>> sitesByGene = groupSitesByGene(sites);
        List<SiteResult> siteResults = testSites(
                samples,
                sitesByGene,
                control,
                caseCondition,
                minTotal,
                geneTotals
        );
        applySiteFdr(siteResults);

        Results results = new Results();
        results.sites = siteResults;

        Map<String, List<SiteResult>> siteResultsByGene =
                groupSiteResultsByGene(siteResults);
        for (Map.Entry<String, List<SiteCluster>> entry : sitesByGene.entrySet()) {
            results.genes.add(summarizeGene(
                    entry.getKey(),
                    entry.getValue(),
                    samples,
                    siteResultsByGene.get(entry.getKey()),
                    control,
                    caseCondition,
                    effect,
                    fdr,
                    wpiThreshold
            ));
        }

        applyGeneFdr(results.genes, effect, fdr, wpiThreshold);
        results.genes.sort((left, right) -> {
            int byKs = Double.compare(left.ksFdr, right.ksFdr);
            return byKs != 0
                    ? byKs
                    : Double.compare(left.minFdr, right.minFdr);
        });
        return results;
    }

    private static List<SiteResult> testSites(
            List<Sample> samples,
            Map<String, List<SiteCluster>> sitesByGene,
            String control,
            String caseCondition,
            int minTotal,
            Map<String, Map<String, Integer>> externalGeneTotals
    ) {
        List<SiteResult> results = new ArrayList<SiteResult>();

        for (Map.Entry<String, List<SiteCluster>> entry : sitesByGene.entrySet()) {
            String geneId = entry.getKey();
            List<SiteCluster> geneSites = entry.getValue();
            int[] geneTotals = geneTotals(
                    samples,
                    geneId,
                    geneSites,
                    externalGeneTotals
            );

            for (SiteCluster site : geneSites) {
                if (site.totalCount() < minTotal) {
                    continue;
                }
                SiteResult result = testSite(
                        samples,
                        site,
                        geneTotals,
                        control,
                        caseCondition
                );
                if (result != null) {
                    results.add(result);
                }
            }
        }
        return results;
    }

    private static SiteResult testSite(
            List<Sample> samples,
            SiteCluster site,
            int[] geneTotals,
            String control,
            String caseCondition
    ) {
        int sampleCount = samples.size();
        int[] siteCounts = new int[sampleCount];
        int[] groups = new int[sampleCount];
        java.util.Arrays.fill(groups, -2);

        double controlPauSum = 0.0;
        double casePauSum = 0.0;
        int controlReplicates = 0;
        int caseReplicates = 0;

        for (int index = 0; index < sampleCount; index++) {
            Sample sample = samples.get(index);
            siteCounts[index] = count(site, sample.name);

            if (sample.condition.equals(control)) {
                groups[index] = 0;
                if (geneTotals[index] > 0) {
                    controlReplicates++;
                    controlPauSum += (double) siteCounts[index] / geneTotals[index];
                }
            } else if (sample.condition.equals(caseCondition)) {
                groups[index] = 1;
                if (geneTotals[index] > 0) {
                    caseReplicates++;
                    casePauSum += (double) siteCounts[index] / geneTotals[index];
                }
            }
        }

        if (controlReplicates == 0 || caseReplicates == 0) {
            return null;
        }

        SiteResult result = new SiteResult();
        result.site = site;
        result.controlPau = controlPauSum / controlReplicates;
        result.casePau = casePauSum / caseReplicates;
        result.delta = result.casePau - result.controlPau;

        if (controlReplicates >= 2 && caseReplicates >= 2) {
            BetaBinomial.TestResult test = BetaBinomial.test(
                    siteCounts,
                    geneTotals,
                    groups
            );
            result.p = test.pvalue;
            result.rho = test.rho;
            result.method = "beta-binomial-LRT";
        } else {
            int[] table = pooledTwoByTwo(siteCounts, geneTotals, groups);
            result.p = FisherExact.twoSided(
                    table[0],
                    table[1],
                    table[2],
                    table[3]
            );
            result.rho = Double.NaN;
            result.method = "Fisher-exact";
        }
        return result;
    }

    private static int[] geneTotals(
            List<Sample> samples,
            String geneId,
            List<SiteCluster> geneSites,
            Map<String, Map<String, Integer>> externalGeneTotals
    ) {
        int[] totals = new int[samples.size()];
        Map<String, Integer> provided = externalGeneTotals == null
                ? null
                : externalGeneTotals.get(geneId);

        for (int index = 0; index < samples.size(); index++) {
            Sample sample = samples.get(index);
            Integer external = provided == null ? null : provided.get(sample.name);
            if (external != null) {
                totals[index] = external;
                continue;
            }
            for (SiteCluster site : geneSites) {
                totals[index] += count(site, sample.name);
            }
        }
        return totals;
    }

    private static int[] pooledTwoByTwo(
            int[] siteCounts,
            int[] geneTotals,
            int[] groups
    ) {
        int controlSite = 0;
        int controlOther = 0;
        int caseSite = 0;
        int caseOther = 0;

        for (int index = 0; index < groups.length; index++) {
            if (groups[index] == 0) {
                controlSite += siteCounts[index];
                controlOther += Math.max(0, geneTotals[index] - siteCounts[index]);
            } else if (groups[index] == 1) {
                caseSite += siteCounts[index];
                caseOther += Math.max(0, geneTotals[index] - siteCounts[index]);
            }
        }
        return new int[]{controlSite, controlOther, caseSite, caseOther};
    }

    private static void applySiteFdr(List<SiteResult> results) {
        double[] pValues = new double[results.size()];
        for (int index = 0; index < results.size(); index++) {
            pValues[index] = results.get(index).p;
        }
        double[] adjusted = MultipleTesting.bh(pValues);
        for (int index = 0; index < results.size(); index++) {
            results.get(index).q = adjusted[index];
        }
    }

    private static GeneResult summarizeGene(
            String geneId,
            List<SiteCluster> geneSites,
            List<Sample> samples,
            List<SiteResult> testedSites,
            String control,
            String caseCondition,
            double effect,
            double fdr,
            double wpiThreshold
    ) {
        geneSites.sort(Comparator.comparingInt(site -> site.siteRank));

        GeneResult result = new GeneResult();
        result.geneId = geneId;
        result.geneName = geneSites.get(0).geneName;
        result.biotype = geneSites.get(0).biotype;
        result.nSites = geneSites.size();
        result.apaType = classifyApaType(geneSites);

        double[] controlMeanPau = new double[geneSites.size()];
        double[] caseMeanPau = new double[geneSites.size()];
        Map<Integer, Integer> controlPositions = new TreeMap<Integer, Integer>();
        Map<Integer, Integer> casePositions = new TreeMap<Integer, Integer>();

        double controlWpiSum = 0.0;
        double caseWpiSum = 0.0;
        int controlReplicates = 0;
        int caseReplicates = 0;

        for (Sample sample : samples) {
            int total = totalCount(geneSites, sample.name);
            if (total == 0) {
                continue;
            }

            double sampleWpi = 0.0;
            for (int siteIndex = 0; siteIndex < geneSites.size(); siteIndex++) {
                SiteCluster site = geneSites.get(siteIndex);
                int siteCount = count(site, sample.name);
                double pau = (double) siteCount / total;
                sampleWpi += pau * site.normalizedPosition;

                if (sample.condition.equals(control)) {
                    controlMeanPau[siteIndex] += pau;
                    addPositionCount(
                            controlPositions,
                            site.representativePos0,
                            siteCount
                    );
                } else if (sample.condition.equals(caseCondition)) {
                    caseMeanPau[siteIndex] += pau;
                    addPositionCount(
                            casePositions,
                            site.representativePos0,
                            siteCount
                    );
                }
            }

            if (sample.condition.equals(control)) {
                controlWpiSum += sampleWpi;
                controlReplicates++;
            } else if (sample.condition.equals(caseCondition)) {
                caseWpiSum += sampleWpi;
                caseReplicates++;
            }
        }

        result.controlWpi = mean(controlWpiSum, controlReplicates);
        result.caseWpi = mean(caseWpiSum, caseReplicates);
        averageInPlace(controlMeanPau, controlReplicates);
        averageInPlace(caseMeanPau, caseReplicates);
        result.deltaWpi = result.caseWpi - result.controlWpi;
        result.apd = totalVariation(controlMeanPau, caseMeanPau);

        TieAdjustedKs.Result ks = TieAdjustedKs.test(
                controlPositions,
                casePositions
        );
        boolean minusStrand = geneSites.get(0).strand == Strand.MINUS;
        result.ksD = ks.d;
        result.ksSignedD = minusStrand ? -ks.signedD : ks.signedD;
        result.ksZ = minusStrand ? -ks.z : ks.z;
        result.ksPvalue = ks.pvalue;
        result.ksMaxDiffAt0 = ks.maxDiffAt0;
        result.ksFdr = 1.0;

        result.minFdr = 1.0;
        boolean significantSiteChange = false;
        if (testedSites != null) {
            for (SiteResult siteResult : testedSites) {
                result.minFdr = Math.min(result.minFdr, siteResult.q);
                if (siteResult.q <= fdr && Math.abs(siteResult.delta) >= effect) {
                    significantSiteChange = true;
                }
            }
        }

        result.apaChange = significantSiteChange;
        result.apaDirection = direction(
                significantSiteChange,
                result.deltaWpi,
                wpiThreshold
        );
        return result;
    }

    private static void applyGeneFdr(
            List<GeneResult> genes,
            double effect,
            double fdr,
            double wpiThreshold
    ) {
        double[] pValues = new double[genes.size()];
        for (int index = 0; index < genes.size(); index++) {
            pValues[index] = genes.get(index).ksPvalue;
        }

        double[] adjusted = MultipleTesting.bh(pValues);
        for (int index = 0; index < genes.size(); index++) {
            GeneResult gene = genes.get(index);
            gene.ksFdr = adjusted[index];

            if (!"stable".equals(gene.apaDirection) || gene.ksFdr > fdr) {
                continue;
            }

            if (gene.deltaWpi >= wpiThreshold) {
                gene.apaDirection = "lengthening";
                gene.apaChange = true;
            } else if (gene.deltaWpi <= -wpiThreshold) {
                gene.apaDirection = "shortening";
                gene.apaChange = true;
            } else if (gene.apd >= effect) {
                gene.apaDirection = "no_directional_shift";
                gene.apaChange = true;
            }
        }
    }

    private static String direction(
            boolean significant,
            double deltaWpi,
            double threshold
    ) {
        if (!significant) {
            return "stable";
        }
        if (deltaWpi >= threshold) {
            return "lengthening";
        }
        if (deltaWpi <= -threshold) {
            return "shortening";
        }
        return "no_directional_shift";
    }

    private static Map<String, List<SiteCluster>> groupSitesByGene(
            List<SiteCluster> sites
    ) {
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

    private static Map<String, List<SiteResult>> groupSiteResultsByGene(
            List<SiteResult> results
    ) {
        Map<String, List<SiteResult>> grouped =
                new HashMap<String, List<SiteResult>>();
        for (SiteResult result : results) {
            grouped.computeIfAbsent(
                    result.site.geneId,
                    key -> new ArrayList<SiteResult>()
            ).add(result);
        }
        return grouped;
    }

    private static String classifyApaType(List<SiteCluster> sites) {
        if (sites.size() <= 1) {
            return "none";
        }

        Set<Integer> terminalAnchors = new LinkedHashSet<Integer>();
        boolean allTerminal = true;
        for (SiteCluster site : sites) {
            if (!site.terminalExon) {
                allTerminal = false;
            } else if (site.terminalExonAnchor0 >= 0) {
                terminalAnchors.add(site.terminalExonAnchor0);
            }
        }

        if (allTerminal && terminalAnchors.size() > 1) {
            return "ALE";
        }
        if (allTerminal && terminalAnchors.size() == 1) {
            return "TUTR";
        }
        return "none";
    }

    private static int totalCount(List<SiteCluster> sites, String sample) {
        int total = 0;
        for (SiteCluster site : sites) {
            total += count(site, sample);
        }
        return total;
    }

    private static int count(SiteCluster site, String sample) {
        Integer value = site.counts.get(sample);
        return value == null ? 0 : value;
    }

    private static void addPositionCount(
            Map<Integer, Integer> counts,
            int position,
            int value
    ) {
        if (value <= 0) {
            return;
        }
        Integer previous = counts.get(position);
        counts.put(position, (previous == null ? 0 : previous) + value);
    }

    private static double mean(double sum, int count) {
        return count > 0 ? sum / count : 0.0;
    }

    private static void averageInPlace(double[] values, int count) {
        if (count <= 0) {
            return;
        }
        for (int index = 0; index < values.length; index++) {
            values[index] /= count;
        }
    }

    private static double totalVariation(double[] left, double[] right) {
        double distance = 0.0;
        for (int index = 0; index < left.length; index++) {
            distance += Math.abs(right[index] - left[index]);
        }
        return 0.5 * distance;
    }
}
