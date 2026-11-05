package org.apaura.analysis;

import org.apaura.model.Sample;
import org.apaura.model.ShortReadGroup;
import org.apaura.model.ShortReadQuantification;
import org.apaura.model.ShortReadSite;
import org.apaura.stats.BetaBinomial;
import org.apaura.stats.FisherExact;
import org.apaura.stats.MultipleTesting;
import org.apaura.stats.ModeratedDirichletMultinomial;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/** Differential analysis for coverage-derived short-read PAS usage. */
public final class ShortReadDifferential {

    private static final int DEFAULT_MIN_INFORMATIVE_REPLICATES = 3;

    private ShortReadDifferential() {
        // Utility class.
    }

    public static final class SiteResult {
        public ShortReadGroup group;
        public ShortReadSite site;

        public double controlPau = Double.NaN;
        public double casePau = Double.NaN;
        public double delta = Double.NaN;
        public double controlPauSd = Double.NaN;
        public double casePauSd = Double.NaN;

        public int controlTotalReplicates;
        public int caseTotalReplicates;
        public int controlAnyEvidence;
        public int caseAnyEvidence;
        public int controlInformativeReplicates;
        public int caseInformativeReplicates;
        public int controlEvidence;
        public int caseEvidence;

        public double p = Double.NaN;
        public double q = Double.NaN;
        public double rho = Double.NaN;
        public String method = "not_tested";
        public String status = "insufficient_evidence";
    }

    public static final class GroupResult {
        public ShortReadGroup group;

        public double controlWpi = Double.NaN;
        public double caseWpi = Double.NaN;
        public double deltaWpi = Double.NaN;
        public double controlWpiSd = Double.NaN;
        public double caseWpiSd = Double.NaN;
        /** Total-variation distance between mean control and case cPAU vectors. */
        public double apd = Double.NaN;
        /** Jensen-Shannon divergence between mean control and case cPAU vectors. */
        public double jsd = Double.NaN;
        /**
         * Formal TUTR-event p-value. This is computed from the replicate-aware
         * PAS p-values within the group using a Bonferroni global-null test.
         */
        public double eventPvalue = Double.NaN;

        /**
         * BH-adjusted TUTR-event FDR across all testable TUTR groups.
         */
        public double eventFdr = Double.NaN;
        public double eventPrecision = Double.NaN;
        public String eventMethod = "not_tested";

        /**
         * Diagnostic minimum PAS-level p-value within this TUTR group.
         */
        public double minSitePvalue = Double.NaN;

        /**
         * Diagnostic minimum PAS-level FDR within this TUTR group.
         */
        public double minSiteFdr = Double.NaN;

        /**
         * Backward-compatible aliases used by existing APAURA writers.
         * These now mirror eventPvalue/eventFdr rather than PAS-level minima.
         */
        public double minPvalue = Double.NaN;
        public double minFdr = Double.NaN;

        public int controlTotalReplicates;
        public int caseTotalReplicates;
        public int controlAnyEvidence;
        public int caseAnyEvidence;
        public int controlInformativeReplicates;
        public int caseInformativeReplicates;

        public String apaArchitecture = "TUTR";
        public String observedDirection = "not_applicable";
        public String apaDirection = "not_applicable";
        public String apaEvent = "TUTR:not_applicable";
        public String status = "insufficient_evidence";
        public boolean apaChange = false;
    }

    public static final class Results {
        public final List<SiteResult> sites = new ArrayList<SiteResult>();
        public final List<GroupResult> groups = new ArrayList<GroupResult>();
    }

    /**
     * Backward-compatible entry point. Replicated experiments require at least
     * three informative biological replicates per condition by default.
     */
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
        return run(
                samples,
                bySample,
                control,
                caseCondition,
                minSampleEvidence,
                DEFAULT_MIN_INFORMATIVE_REPLICATES,
                effect,
                fdr,
                wpiThreshold
        );
    }

    /**
     * Replicate-aware differential analysis.
     *
     * <p>For replicated experiments, a PAS is statistically tested only when
     * both conditions have at least {@code minInformativeReplicates} samples
     * with {@code commonEvidenceReads >= minSampleEvidence}. Lack of replicate
     * support is reported explicitly and is not rescued by pooling reads into
     * Fisher's exact test. Fisher's exact test is retained only for genuine
     * 1-vs-1 analyses.</p>
     */
    public static Results run(
            List<Sample> samples,
            Map<String, List<ShortReadQuantification>> bySample,
            String control,
            String caseCondition,
            int minSampleEvidence,
            int minInformativeReplicates,
            double effect,
            double fdr,
            double wpiThreshold
    ) {
        if (minInformativeReplicates < 1) {
            throw new IllegalArgumentException("minInformativeReplicates must be >= 1");
        }

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
                minInformativeReplicates,
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
                minInformativeReplicates,
                effect,
                fdr,
                wpiThreshold,
                results
        );

        applyModeratedDmToGroups(
                samples,
                bySample,
                control,
                caseCondition,
                minSampleEvidence,
                results.groups
        );
        applyGroupFdr(results.groups);
        classifyGroups(results.groups, effect, fdr, wpiThreshold);

        results.sites.sort((left, right) -> compareFinite(left.q, right.q));
        results.groups.sort((left, right) -> compareFinite(left.eventFdr, right.eventFdr));
        return results;
    }

    private static void testSites(
            List<Sample> samples,
            Map<String, List<ShortReadQuantification>> bySample,
            List<ShortReadQuantification> firstSample,
            String control,
            String caseCondition,
            int minSampleEvidence,
            int minInformativeReplicates,
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
                        minSampleEvidence,
                        minInformativeReplicates
                );

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
            int minSampleEvidence,
            int minInformativeReplicates
    ) {
        ShortReadGroup group = bySample.get(samples.get(0).name)
                .get(groupIndex)
                .group;

        int[] siteCounts = new int[samples.size()];
        int[] totalCounts = new int[samples.size()];
        int[] groups = new int[samples.size()];
        Arrays.fill(groups, -2);

        List<Double> controlPau = new ArrayList<Double>();
        List<Double> casePau = new ArrayList<Double>();

        SiteResult result = new SiteResult();
        result.group = group;
        result.site = group.sites.get(siteIndex);

        for (int index = 0; index < samples.size(); index++) {
            Sample sample = samples.get(index);
            ShortReadQuantification quantification =
                    bySample.get(sample.name).get(groupIndex);

            totalCounts[index] = quantification.commonEvidenceReads;
            siteCounts[index] = quantification.pseudoCounts[siteIndex];

            boolean isControl = sample.condition.equals(control);
            boolean isCase = sample.condition.equals(caseCondition);
            if (!isControl && !isCase) {
                continue;
            }

            if (isControl) {
                result.controlTotalReplicates++;
                if (totalCounts[index] > 0) {
                    result.controlAnyEvidence++;
                }
            } else {
                result.caseTotalReplicates++;
                if (totalCounts[index] > 0) {
                    result.caseAnyEvidence++;
                }
            }

            if (totalCounts[index] < minSampleEvidence) {
                continue;
            }

            if (isControl) {
                groups[index] = 0;
                result.controlInformativeReplicates++;
                result.controlEvidence += totalCounts[index];
                controlPau.add(quantification.sitePau[siteIndex]);
            } else {
                groups[index] = 1;
                result.caseInformativeReplicates++;
                result.caseEvidence += totalCounts[index];
                casePau.add(quantification.sitePau[siteIndex]);
            }
        }

        result.controlPau = ApaUsageMetrics.mean(controlPau);
        result.casePau = ApaUsageMetrics.mean(casePau);
        result.controlPauSd = ApaUsageMetrics.sampleSd(controlPau);
        result.casePauSd = ApaUsageMetrics.sampleSd(casePau);
        result.delta = finiteDifference(result.casePau, result.controlPau);

        if (result.controlInformativeReplicates == 0
                || result.caseInformativeReplicates == 0) {
            result.status = "insufficient_evidence";
            return result;
        }

        boolean unreplicated = result.controlTotalReplicates == 1
                && result.caseTotalReplicates == 1;

        if (unreplicated) {
            int[] table = pooledTwoByTwo(siteCounts, totalCounts, groups);
            result.p = FisherExact.twoSided(
                    table[0],
                    table[1],
                    table[2],
                    table[3]
            );
            result.method = "coverage-Fisher-exact-unreplicated";
            result.status = "ok";
            return result;
        }

        if (result.controlInformativeReplicates < minInformativeReplicates
                || result.caseInformativeReplicates < minInformativeReplicates) {
            result.status = "insufficient_replicate_support";
            return result;
        }

        BetaBinomial.TestResult test = BetaBinomial.test(
                siteCounts,
                totalCounts,
                groups
        );
        result.p = test.pvalue;
        result.rho = test.rho;
        result.method = "coverage-beta-binomial-LRT";
        result.status = "ok";
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

    /** Apply BH only to testable finite p-values. Untestable rows remain NA. */
    private static void applySiteFdr(List<SiteResult> results) {
        List<Integer> testedIndices = new ArrayList<Integer>();
        List<Double> testedPValues = new ArrayList<Double>();

        for (int index = 0; index < results.size(); index++) {
            double p = results.get(index).p;
            if (Double.isFinite(p)) {
                testedIndices.add(index);
                testedPValues.add(p);
            }
        }

        if (testedIndices.isEmpty()) {
            return;
        }

        double[] pValues = new double[testedPValues.size()];
        for (int index = 0; index < testedPValues.size(); index++) {
            pValues[index] = testedPValues.get(index);
        }

        double[] adjusted = MultipleTesting.bh(pValues);
        for (int index = 0; index < testedIndices.size(); index++) {
            results.get(testedIndices.get(index)).q = adjusted[index];
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
            int minInformativeReplicates,
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
                    minInformativeReplicates,
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
            int minInformativeReplicates,
            double effect,
            double fdr,
            double wpiThreshold
    ) {
        int siteCount = group.sites.size();
        double[] controlMeanPau = new double[siteCount];
        double[] caseMeanPau = new double[siteCount];

        List<Double> controlWpi = new ArrayList<Double>();
        List<Double> caseWpi = new ArrayList<Double>();

        GroupResult result = new GroupResult();
        result.group = group;

        for (Sample sample : samples) {
            ShortReadQuantification quantification =
                    bySample.get(sample.name).get(groupIndex);

            boolean isControl = sample.condition.equals(control);
            boolean isCase = sample.condition.equals(caseCondition);
            if (!isControl && !isCase) {
                continue;
            }

            if (isControl) {
                result.controlTotalReplicates++;
                if (quantification.commonEvidenceReads > 0) {
                    result.controlAnyEvidence++;
                }
            } else {
                result.caseTotalReplicates++;
                if (quantification.commonEvidenceReads > 0) {
                    result.caseAnyEvidence++;
                }
            }

            if (quantification.commonEvidenceReads < minSampleEvidence) {
                continue;
            }

            double sampleWpi = weightedPositionIndex(group, quantification.sitePau);
            if (isControl) {
                result.controlInformativeReplicates++;
                controlWpi.add(sampleWpi);
                ApaUsageMetrics.addInPlace(controlMeanPau, quantification.sitePau);
            } else {
                result.caseInformativeReplicates++;
                caseWpi.add(sampleWpi);
                ApaUsageMetrics.addInPlace(caseMeanPau, quantification.sitePau);
            }
        }

        if (!controlWpi.isEmpty()) {
            result.controlWpi = ApaUsageMetrics.mean(controlWpi);
            result.controlWpiSd = ApaUsageMetrics.sampleSd(controlWpi);
            ApaUsageMetrics.averageInPlace(controlMeanPau, controlWpi.size());
        }
        if (!caseWpi.isEmpty()) {
            result.caseWpi = ApaUsageMetrics.mean(caseWpi);
            result.caseWpiSd = ApaUsageMetrics.sampleSd(caseWpi);
            ApaUsageMetrics.averageInPlace(caseMeanPau, caseWpi.size());
        }

        result.deltaWpi = finiteDifference(result.caseWpi, result.controlWpi);
        if (!controlWpi.isEmpty() && !caseWpi.isEmpty()) {
            result.apd = ApaUsageMetrics.totalVariationDistance(
                    controlMeanPau,
                    caseMeanPau
            );
            result.jsd = ApaUsageMetrics.jensenShannonDivergence(
                    controlMeanPau,
                    caseMeanPau
            );
            result.observedDirection = ApaUsageMetrics.observedDirection(
                    result.deltaWpi,
                    wpiThreshold
            );
        }

        if (result.controlInformativeReplicates == 0
                || result.caseInformativeReplicates == 0) {
            result.status = "insufficient_evidence";
            return result;
        }

        boolean unreplicated = result.controlTotalReplicates == 1
                && result.caseTotalReplicates == 1;
        if (!unreplicated
                && (result.controlInformativeReplicates < minInformativeReplicates
                || result.caseInformativeReplicates < minInformativeReplicates)) {
            result.status = "insufficient_replicate_support";
            return result;
        }

        result.status = "ok";
        int testedSiteCount = 0;

        if (siteResults != null) {
            for (SiteResult site : siteResults) {
                if ("ok".equals(site.status) && Double.isFinite(site.p)) {
                    testedSiteCount++;

                    if (!Double.isFinite(result.minSitePvalue)
                            || site.p < result.minSitePvalue) {
                        result.minSitePvalue = site.p;
                    }
                }

                if (Double.isFinite(site.q)
                        && (!Double.isFinite(result.minSiteFdr)
                        || site.q < result.minSiteFdr)) {
                    result.minSiteFdr = site.q;
                }
            }
        }

        if (unreplicated) {
            if (testedSiteCount == 0 || !Double.isFinite(result.minSitePvalue)) {
                result.status = "insufficient_evidence";
                return result;
            }

            result.eventPvalue = Math.min(
                    1.0,
                    result.minSitePvalue * testedSiteCount
            );
            result.eventMethod = "bonferroni-min-site-unreplicated";
            result.minPvalue = result.eventPvalue;
        }

        return result;
    }

    private static void applyModeratedDmToGroups(
            List<Sample> samples,
            Map<String, List<ShortReadQuantification>> bySample,
            String control,
            String caseCondition,
            int minSampleEvidence,
            List<GroupResult> groups
    ) {
        List<ModeratedDirichletMultinomial.Event> events =
                new ArrayList<ModeratedDirichletMultinomial.Event>();
        Map<String, GroupResult> byId =
                new HashMap<String, GroupResult>();

        for (int groupIndex = 0; groupIndex < groups.size(); groupIndex++) {
            GroupResult result = groups.get(groupIndex);
            if (!"ok".equals(result.status)) continue;

            boolean unreplicated = result.controlTotalReplicates == 1
                    && result.caseTotalReplicates == 1;
            if (unreplicated) continue;

            List<int[]> rows = new ArrayList<int[]>();
            List<Integer> labels = new ArrayList<Integer>();

            for (Sample sample : samples) {
                boolean isControl = sample.condition.equals(control);
                boolean isCase = sample.condition.equals(caseCondition);
                if (!isControl && !isCase) continue;

                ShortReadQuantification q =
                        bySample.get(sample.name).get(groupIndex);
                if (q.commonEvidenceReads < minSampleEvidence) continue;

                rows.add(Arrays.copyOf(q.pseudoCounts, q.pseudoCounts.length));
                labels.add(isControl ? 0 : 1);
            }

            if (rows.isEmpty()) continue;

            int[][] counts = new int[rows.size()][];
            int[] groupLabels = new int[rows.size()];
            for (int i = 0; i < rows.size(); i++) {
                counts[i] = rows.get(i);
                groupLabels[i] = labels.get(i);
            }

            result.eventPvalue = Double.NaN;
            result.eventFdr = Double.NaN;
            result.minPvalue = Double.NaN;
            result.minFdr = Double.NaN;
            result.eventMethod = "moderated-DM-LRT";

            events.add(new ModeratedDirichletMultinomial.Event(
                    result.group.id,
                    counts,
                    groupLabels
            ));
            byId.put(result.group.id, result);
        }

        ModeratedDirichletMultinomial.BatchResult batch =
                ModeratedDirichletMultinomial.test(events);

        for (ModeratedDirichletMultinomial.Result dm : batch.results) {
            GroupResult result = byId.get(dm.id);
            if (result == null) continue;

            if (dm.converged && Double.isFinite(dm.pvalue)) {
                result.eventPvalue = dm.pvalue;
                result.minPvalue = dm.pvalue;
                result.eventPrecision = dm.precision;
            } else {
                result.status = "dm_fit_failed";
                result.eventMethod = "moderated-DM-LRT-failed";
            }
        }
    }

    /**
     * Apply BH correction across testable TUTR-event p-values.
     *
     * Untestable groups remain NA and are excluded from the multiple-testing
     * denominator, matching the behavior of PAS-level FDR correction.
     */
    private static void applyGroupFdr(List<GroupResult> groups) {
        List<Integer> testedIndices = new ArrayList<Integer>();
        List<Double> testedPValues = new ArrayList<Double>();

        for (int index = 0; index < groups.size(); index++) {
            double p = groups.get(index).eventPvalue;

            if (Double.isFinite(p)) {
                testedIndices.add(index);
                testedPValues.add(p);
            }
        }

        if (testedIndices.isEmpty()) {
            return;
        }

        double[] pValues = new double[testedPValues.size()];

        for (int index = 0; index < testedPValues.size(); index++) {
            pValues[index] = testedPValues.get(index);
        }

        double[] adjusted = MultipleTesting.bh(pValues);

        for (int index = 0; index < testedIndices.size(); index++) {
            GroupResult result = groups.get(testedIndices.get(index));

            result.eventFdr = adjusted[index];

            // Preserve the historical output field while correcting its meaning.
            result.minFdr = result.eventFdr;
        }
    }

    /**
     * Classify biological TUTR remodeling only after event-level FDR has been
     * calculated.
     */
    private static void classifyGroups(
            List<GroupResult> groups,
            double effect,
            double fdr,
            double wpiThreshold
    ) {
        for (GroupResult result : groups) {
            if (!"ok".equals(result.status)) {
                continue;
            }

            boolean statisticallySupported =
                    Double.isFinite(result.eventFdr)
                    && result.eventFdr <= fdr;

            boolean observedChange =
                    Double.isFinite(result.apd)
                    && result.apd >= effect;

            result.apaChange =
                    statisticallySupported && observedChange;

            result.apaDirection = ApaUsageMetrics.supportedDirection(
                    observedChange,
                    result.deltaWpi,
                    wpiThreshold
            );

            result.apaEvent =
                    result.apaArchitecture + ":" + result.apaDirection;
        }
    }

    private static double weightedPositionIndex(
            ShortReadGroup group,
            double[] pau
    ) {
        double[] positions = new double[group.sites.size()];
        for (int index = 0; index < group.sites.size(); index++) {
            positions[index] = group.sites.get(index).normalizedPosition;
        }
        return ApaUsageMetrics.weightedPositionIndex(pau, positions);
    }

    private static double finiteDifference(double right, double left) {
        if (!Double.isFinite(right) || !Double.isFinite(left)) {
            return Double.NaN;
        }
        return right - left;
    }

    private static int compareFinite(double left, double right) {
        boolean leftFinite = Double.isFinite(left);
        boolean rightFinite = Double.isFinite(right);
        if (leftFinite && rightFinite) {
            return Double.compare(left, right);
        }
        if (leftFinite) {
            return -1;
        }
        if (rightFinite) {
            return 1;
        }
        return 0;
    }
}
