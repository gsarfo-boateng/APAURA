package org.apaura.analysis;

import org.apaura.model.Sample;
import org.apaura.model.ShortReadGroup;
import org.apaura.model.ShortReadQuantification;
import org.apaura.model.ShortReadSite;
import org.apaura.stats.BetaBinomial;
import org.apaura.stats.FisherExact;
import org.apaura.stats.MultipleTesting;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/** Differential analysis for coverage-derived short-read PAS usage. */
public final class ShortReadDifferential {

    private ShortReadDifferential() {
        // Utility class.
    }

    public static final class SiteResult {
        public ShortReadGroup group;
        public ShortReadSite site;
        public double controlPau;
        public double casePau;
        public double delta;
        public double p;
        public double q;
        public double rho;
        public String method;
        public int controlEvidence;
        public int caseEvidence;
    }

    public static final class GroupResult {
        public ShortReadGroup group;
        public double controlWpi;
        public double caseWpi;
        public double deltaWpi;
        public double apd;
        public double minFdr;
        public String apaArchitecture = "TUTR";
        public String apaDirection = "stable";
        public String apaEvent = "TUTR:stable";
        public String status;
        public boolean apaChange = false;
    }

    public static final class Results {
        public final List<SiteResult> sites = new ArrayList<SiteResult>();
        public final List<GroupResult> groups = new ArrayList<GroupResult>();
    }

    public static Results run(
            List<Sample> samples,
            Map<String, List<ShortReadQuantification>> bySample,
            String control,
            String caseCondition,
            int minSampleEvidence,
            double effect,
            double fdr,
            double wpiThreshold
    ) {
        Results results = new Results();
        if (samples.isEmpty()) {
            return results;
        }

        List<ShortReadQuantification> firstSample = bySample.get(samples.get(0).name);
        if (firstSample == null) {
            return results;
        }

        Map<String, List<SiteResult>> siteResultsByGroup =
                new HashMap<String, List<SiteResult>>();
        testSites(
                samples,
                bySample,
                firstSample,
                control,
                caseCondition,
                minSampleEvidence,
                results,
                siteResultsByGroup
        );
        applySiteFdr(results.sites);

        summarizeGroups(
                samples,
                bySample,
                firstSample,
                siteResultsByGroup,
                control,
                caseCondition,
                minSampleEvidence,
                effect,
                fdr,
                wpiThreshold,
                results
        );

        results.sites.sort((left, right) -> Double.compare(left.q, right.q));
        results.groups.sort((left, right) ->
                Double.compare(left.minFdr, right.minFdr));
        return results;
    }

    private static void testSites(
            List<Sample> samples,
            Map<String, List<ShortReadQuantification>> bySample,
            List<ShortReadQuantification> firstSample,
            String control,
            String caseCondition,
            int minSampleEvidence,
            Results results,
            Map<String, List<SiteResult>> siteResultsByGroup
    ) {
        for (int groupIndex = 0; groupIndex < firstSample.size(); groupIndex++) {
            ShortReadGroup group = firstSample.get(groupIndex).group;
            for (int siteIndex = 0; siteIndex < group.sites.size(); siteIndex++) {
                SiteResult result = testSite(
                        samples,
                        bySample,
                        groupIndex,
                        siteIndex,
                        control,
                        caseCondition,
                        minSampleEvidence
                );
                if (result == null) {
                    continue;
                }
                results.sites.add(result);
                siteResultsByGroup.computeIfAbsent(
                        group.id,
                        key -> new ArrayList<SiteResult>()
                ).add(result);
            }
        }
    }

    private static SiteResult testSite(
            List<Sample> samples,
            Map<String, List<ShortReadQuantification>> bySample,
            int groupIndex,
            int siteIndex,
            String control,
            String caseCondition,
            int minSampleEvidence
    ) {
        ShortReadGroup group = bySample.get(samples.get(0).name)
                .get(groupIndex)
                .group;

        int[] siteCounts = new int[samples.size()];
        int[] totalCounts = new int[samples.size()];
        int[] groups = new int[samples.size()];
        Arrays.fill(groups, -2);

        double controlPauSum = 0.0;
        double casePauSum = 0.0;
        int controlReplicates = 0;
        int caseReplicates = 0;
        int controlEvidence = 0;
        int caseEvidence = 0;

        for (int index = 0; index < samples.size(); index++) {
            Sample sample = samples.get(index);
            ShortReadQuantification quantification =
                    bySample.get(sample.name).get(groupIndex);

            totalCounts[index] = quantification.commonEvidenceReads;
            siteCounts[index] = quantification.pseudoCounts[siteIndex];
            if (totalCounts[index] < minSampleEvidence) {
                continue;
            }

            if (sample.condition.equals(control)) {
                groups[index] = 0;
                controlPauSum += quantification.sitePau[siteIndex];
                controlReplicates++;
                controlEvidence += totalCounts[index];
            } else if (sample.condition.equals(caseCondition)) {
                groups[index] = 1;
                casePauSum += quantification.sitePau[siteIndex];
                caseReplicates++;
                caseEvidence += totalCounts[index];
            }
        }

        if (controlReplicates == 0 || caseReplicates == 0) {
            return null;
        }

        SiteResult result = new SiteResult();
        result.group = group;
        result.site = group.sites.get(siteIndex);
        result.controlPau = controlPauSum / controlReplicates;
        result.casePau = casePauSum / caseReplicates;
        result.delta = result.casePau - result.controlPau;
        result.controlEvidence = controlEvidence;
        result.caseEvidence = caseEvidence;

        if (controlReplicates >= 2 && caseReplicates >= 2) {
            BetaBinomial.TestResult test = BetaBinomial.test(
                    siteCounts,
                    totalCounts,
                    groups
            );
            result.p = test.pvalue;
            result.rho = test.rho;
            result.method = "coverage-beta-binomial-LRT";
        } else {
            int[] table = pooledTwoByTwo(siteCounts, totalCounts, groups);
            result.p = FisherExact.twoSided(
                    table[0],
                    table[1],
                    table[2],
                    table[3]
            );
            result.rho = Double.NaN;
            result.method = "coverage-Fisher-exact";
        }
        return result;
    }

    private static int[] pooledTwoByTwo(
            int[] siteCounts,
            int[] totalCounts,
            int[] groups
    ) {
        int controlSite = 0;
        int controlOther = 0;
        int caseSite = 0;
        int caseOther = 0;

        for (int index = 0; index < groups.length; index++) {
            if (groups[index] == 0) {
                controlSite += siteCounts[index];
                controlOther += totalCounts[index] - siteCounts[index];
            } else if (groups[index] == 1) {
                caseSite += siteCounts[index];
                caseOther += totalCounts[index] - siteCounts[index];
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

    private static void summarizeGroups(
            List<Sample> samples,
            Map<String, List<ShortReadQuantification>> bySample,
            List<ShortReadQuantification> firstSample,
            Map<String, List<SiteResult>> siteResultsByGroup,
            String control,
            String caseCondition,
            int minSampleEvidence,
            double effect,
            double fdr,
            double wpiThreshold,
            Results results
    ) {
        for (int groupIndex = 0; groupIndex < firstSample.size(); groupIndex++) {
            ShortReadGroup group = firstSample.get(groupIndex).group;
            GroupResult result = summarizeGroup(
                    samples,
                    bySample,
                    group,
                    groupIndex,
                    siteResultsByGroup.get(group.id),
                    control,
                    caseCondition,
                    minSampleEvidence,
                    effect,
                    fdr,
                    wpiThreshold
            );
            results.groups.add(result);
        }
    }

    private static GroupResult summarizeGroup(
            List<Sample> samples,
            Map<String, List<ShortReadQuantification>> bySample,
            ShortReadGroup group,
            int groupIndex,
            List<SiteResult> siteResults,
            String control,
            String caseCondition,
            int minSampleEvidence,
            double effect,
            double fdr,
            double wpiThreshold
    ) {
        int siteCount = group.sites.size();
        double[] controlMeanPau = new double[siteCount];
        double[] caseMeanPau = new double[siteCount];
        double controlWpiSum = 0.0;
        double caseWpiSum = 0.0;
        int controlReplicates = 0;
        int caseReplicates = 0;

        for (Sample sample : samples) {
            ShortReadQuantification quantification =
                    bySample.get(sample.name).get(groupIndex);
            if (quantification.commonEvidenceReads < minSampleEvidence) {
                continue;
            }

            double sampleWpi = weightedPositionIndex(group, quantification.sitePau);
            if (sample.condition.equals(control)) {
                controlWpiSum += sampleWpi;
                controlReplicates++;
                addInPlace(controlMeanPau, quantification.sitePau);
            } else if (sample.condition.equals(caseCondition)) {
                caseWpiSum += sampleWpi;
                caseReplicates++;
                addInPlace(caseMeanPau, quantification.sitePau);
            }
        }

        GroupResult result = new GroupResult();
        result.group = group;
        result.minFdr = 1.0;

        if (controlReplicates == 0 || caseReplicates == 0) {
            result.status = "insufficient_evidence";
            result.apaDirection = "not_applicable";
            result.apaEvent = "TUTR:not_applicable";
            return result;
        }

        result.controlWpi = controlWpiSum / controlReplicates;
        result.caseWpi = caseWpiSum / caseReplicates;
        averageInPlace(controlMeanPau, controlReplicates);
        averageInPlace(caseMeanPau, caseReplicates);
        result.deltaWpi = result.caseWpi - result.controlWpi;
        result.apd = totalVariation(controlMeanPau, caseMeanPau);
        result.status = "ok";

        boolean significant = false;
        if (siteResults != null) {
            for (SiteResult site : siteResults) {
                result.minFdr = Math.min(result.minFdr, site.q);
                if (site.q <= fdr && Math.abs(site.delta) >= effect) {
                    significant = true;
                }
            }
        }

        result.apaChange = significant;
        result.apaDirection = direction(
                significant,
                result.deltaWpi,
                wpiThreshold
        );
        result.apaEvent = result.apaArchitecture + ":" + result.apaDirection;
        return result;
    }

    private static double weightedPositionIndex(
            ShortReadGroup group,
            double[] pau
    ) {
        double wpi = 0.0;
        for (int index = 0; index < group.sites.size(); index++) {
            wpi += pau[index] * group.sites.get(index).normalizedPosition;
        }
        return wpi;
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

    private static void addInPlace(double[] target, double[] values) {
        for (int index = 0; index < target.length; index++) {
            target[index] += values[index];
        }
    }

    private static void averageInPlace(double[] values, int count) {
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
