package org.apaura.analysis;

import org.apaura.model.CoverageSegment;
import org.apaura.model.Sample;
import org.apaura.model.ShortReadGroup;
import org.apaura.model.ShortReadQuantification;
import org.apaura.model.ShortReadSite;
import org.apaura.model.Strand;
import org.apaura.stats.MultipleTesting;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

class ShortReadDifferentialTest {

    @Test
    void doesNotRescueSparseReplicatesWithPooledFisher() {
        ShortReadGroup group = group();
        List<Sample> samples = new ArrayList<Sample>();
        Map<String, List<ShortReadQuantification>> bySample =
                new LinkedHashMap<String, List<ShortReadQuantification>>();

        for (int i = 1; i <= 10; i++) {
            String name = "WT" + i;
            samples.add(new Sample(name, "control", name + ".bam"));
            bySample.put(name, singleton(quant(group, 0, 0, 0, 0)));
        }
        for (int i = 1; i <= 10; i++) {
            String name = "MUT" + i;
            samples.add(new Sample(name, "case", name + ".bam"));
            bySample.put(name, singleton(quant(group, 0, 0, 0, 0)));
        }

        // GRM1-like support pattern: one informative WT and two informative MUT.
        bySample.put("WT8", singleton(quant(group, 24, 0, 2, 26)));
        bySample.put("WT10", singleton(quant(group, 1, 0, 0, 1)));
        bySample.put("MUT2", singleton(quant(group, 0, 0, 53, 53)));
        bySample.put("MUT4", singleton(quant(group, 2, 0, 0, 2)));
        bySample.put("MUT6", singleton(quant(group, 0, 0, 40, 40)));

        ShortReadDifferential.Results results = ShortReadDifferential.run(
                samples,
                bySample,
                "control",
                "case",
                10,
                3,
                0.10,
                0.05,
                0.10
        );

        ShortReadDifferential.SiteResult distal = results.sites.stream()
                .filter(x -> x.site.id.endsWith("PAS3"))
                .findFirst()
                .orElseThrow(AssertionError::new);

        assertEquals(10, distal.controlTotalReplicates);
        assertEquals(10, distal.caseTotalReplicates);
        assertEquals(2, distal.controlAnyEvidence);
        assertEquals(3, distal.caseAnyEvidence);
        assertEquals(1, distal.controlInformativeReplicates);
        assertEquals(2, distal.caseInformativeReplicates);
        assertEquals("insufficient_replicate_support", distal.status);
        assertEquals("not_tested", distal.method);
        assertTrue(Double.isNaN(distal.p));
        assertTrue(Double.isNaN(distal.q));

        ShortReadDifferential.GroupResult summary = results.groups.get(0);
        assertEquals("insufficient_replicate_support", summary.status);
        assertEquals("lengthening", summary.observedDirection);
        assertEquals("not_applicable", summary.apaDirection);
        assertFalse(summary.apaChange);
        assertTrue(Double.isNaN(summary.eventPvalue));
        assertTrue(Double.isNaN(summary.eventFdr));
        assertTrue(Double.isNaN(summary.minPvalue));
        assertTrue(Double.isNaN(summary.minFdr));
    }

    @Test
    void usesBetaBinomialWhenBothConditionsHaveEnoughReplicates() {
        ShortReadGroup group = group();
        List<Sample> samples = new ArrayList<Sample>();
        Map<String, List<ShortReadQuantification>> bySample =
                new LinkedHashMap<String, List<ShortReadQuantification>>();

        for (int i = 1; i <= 3; i++) {
            String name = "WT" + i;
            samples.add(new Sample(name, "control", name + ".bam"));
            bySample.put(name, singleton(quant(group, 18, 0, 2, 20)));
        }
        for (int i = 1; i <= 3; i++) {
            String name = "MUT" + i;
            samples.add(new Sample(name, "case", name + ".bam"));
            bySample.put(name, singleton(quant(group, 2, 0, 18, 20)));
        }

        ShortReadDifferential.Results results = ShortReadDifferential.run(
                samples,
                bySample,
                "control",
                "case",
                10,
                3,
                0.10,
                0.05,
                0.10
        );

        ShortReadDifferential.SiteResult distal = results.sites.stream()
                .filter(x -> x.site.id.endsWith("PAS3"))
                .findFirst()
                .orElseThrow(AssertionError::new);

        assertEquals(3, distal.controlInformativeReplicates);
        assertEquals(3, distal.caseInformativeReplicates);
        assertEquals("ok", distal.status);
        assertEquals("coverage-beta-binomial-LRT", distal.method);
        assertTrue(Double.isFinite(distal.p));
        assertTrue(Double.isFinite(distal.q));
    }


    @Test
    void detectsMultiPasRedistributionFromVectorMagnitude() {
        ShortReadGroup group = fourSiteGroup();
        List<Sample> samples = Arrays.asList(
                new Sample("C1", "control", "C1.bam"),
                new Sample("T1", "case", "T1.bam")
        );

        Map<String, List<ShortReadQuantification>> bySample =
                new LinkedHashMap<String, List<ShortReadQuantification>>();
        bySample.put(
                "C1",
                singleton(quant4(group, 250, 250, 250, 250))
        );
        bySample.put(
                "T1",
                singleton(quant4(group, 190, 310, 190, 310))
        );

        ShortReadDifferential.Results results = ShortReadDifferential.run(
                samples,
                bySample,
                "control",
                "case",
                10,
                3,
                0.10,
                0.05,
                0.10
        );

        for (ShortReadDifferential.SiteResult site : results.sites) {
            assertTrue(Math.abs(site.delta) < 0.10);
            assertEquals("ok", site.status);
        }

        ShortReadDifferential.GroupResult summary = results.groups.get(0);
        assertEquals("ok", summary.status);
        assertEquals(0.12, summary.apd, 1.0E-12);
        assertTrue(summary.jsd > 0.0);
        assertEquals(0.04, summary.deltaWpi, 1.0E-12);
        assertEquals("redistribution", summary.observedDirection);
        assertTrue(summary.apaChange);
        assertEquals("redistribution", summary.apaDirection);
    }


    @Test
    void computesBonferroniTutrEventPvalueAndKeepsSiteDiagnostics() {
        ShortReadGroup group = group();

        List<Sample> samples = Arrays.asList(
                new Sample("C1", "control", "C1.bam"),
                new Sample("T1", "case", "T1.bam")
        );

        Map<String, List<ShortReadQuantification>> bySample =
                new LinkedHashMap<String, List<ShortReadQuantification>>();

        bySample.put(
                "C1",
                singleton(quant(group, 80, 10, 10, 100))
        );

        bySample.put(
                "T1",
                singleton(quant(group, 10, 10, 80, 100))
        );

        ShortReadDifferential.Results results =
                ShortReadDifferential.run(
                        samples,
                        bySample,
                        "control",
                        "case",
                        10,
                        3,
                        0.10,
                        0.05,
                        0.10
                );

        ShortReadDifferential.GroupResult summary =
                groupResult(results, group.id);

        int testedSites = 0;
        double minSiteP = Double.POSITIVE_INFINITY;
        double minSiteQ = Double.POSITIVE_INFINITY;

        for (ShortReadDifferential.SiteResult site : results.sites) {
            if (!group.id.equals(site.group.id)) {
                continue;
            }

            if ("ok".equals(site.status)
                    && Double.isFinite(site.p)) {
                testedSites++;
                minSiteP = Math.min(minSiteP, site.p);
            }

            if (Double.isFinite(site.q)) {
                minSiteQ = Math.min(minSiteQ, site.q);
            }
        }

        assertEquals(3, testedSites);

        double expectedEventP =
                Math.min(1.0, minSiteP * testedSites);

        /*
         * There is only one TUTR event in this test, so BH(event p)
         * must equal the event p-value itself.
         */
        assertEquals(
                expectedEventP,
                summary.eventPvalue,
                1.0E-12
        );

        assertEquals(
                expectedEventP,
                summary.eventFdr,
                1.0E-12
        );

        assertEquals(
                minSiteP,
                summary.minSitePvalue,
                1.0E-12
        );

        assertEquals(
                minSiteQ,
                summary.minSiteFdr,
                1.0E-12
        );

        /*
         * Historical aliases must now represent event-level inference,
         * not the minimum PAS-level statistics.
         */
        assertEquals(
                summary.eventPvalue,
                summary.minPvalue,
                1.0E-12
        );

        assertEquals(
                summary.eventFdr,
                summary.minFdr,
                1.0E-12
        );

        assertTrue(Double.isFinite(summary.eventPvalue));
        assertTrue(Double.isFinite(summary.eventFdr));
    }


    @Test
    void appliesBhAcrossTutrEventsRatherThanPasTests() {
        ShortReadGroup group1 = group();
        ShortReadGroup group2 = fourSiteGroup();

        List<Sample> samples = Arrays.asList(
                new Sample("C1", "control", "C1.bam"),
                new Sample("T1", "case", "T1.bam")
        );

        Map<String, List<ShortReadQuantification>> bySample =
                new LinkedHashMap<String, List<ShortReadQuantification>>();

        bySample.put(
                "C1",
                Arrays.asList(
                        quant(group1, 80, 10, 10, 100),
                        quant4(group2, 25, 25, 25, 25)
                )
        );

        bySample.put(
                "T1",
                Arrays.asList(
                        quant(group1, 10, 10, 80, 100),
                        quant4(group2, 25, 25, 25, 25)
                )
        );

        ShortReadDifferential.Results results =
                ShortReadDifferential.run(
                        samples,
                        bySample,
                        "control",
                        "case",
                        10,
                        3,
                        0.10,
                        0.05,
                        0.10
                );

        ShortReadDifferential.GroupResult result1 =
                groupResult(results, group1.id);

        ShortReadDifferential.GroupResult result2 =
                groupResult(results, group2.id);

        assertTrue(Double.isFinite(result1.eventPvalue));
        assertTrue(Double.isFinite(result2.eventPvalue));

        double[] eventPvalues = new double[]{
                result1.eventPvalue,
                result2.eventPvalue
        };

        double[] expectedEventFdr =
                MultipleTesting.bh(eventPvalues);

        assertEquals(
                expectedEventFdr[0],
                result1.eventFdr,
                1.0E-12
        );

        assertEquals(
                expectedEventFdr[1],
                result2.eventFdr,
                1.0E-12
        );

        /*
         * The backwards-compatible FDR field must follow the TUTR-event
         * BH correction, not the PAS/site-level BH correction.
         */
        assertEquals(
                result1.eventFdr,
                result1.minFdr,
                1.0E-12
        );

        assertEquals(
                result2.eventFdr,
                result2.minFdr,
                1.0E-12
        );
    }


    private static ShortReadDifferential.GroupResult groupResult(
            ShortReadDifferential.Results results,
            String groupId
    ) {
        for (ShortReadDifferential.GroupResult result : results.groups) {
            if (groupId.equals(result.group.id)) {
                return result;
            }
        }

        throw new AssertionError(
                "Missing TUTR group result: " + groupId
        );
    }


    private static ShortReadGroup group() {
        List<ShortReadSite> sites = Arrays.asList(
                new ShortReadSite("G1|PAS1", 100, 1, 0.0),
                new ShortReadSite("G1|PAS2", 200, 2, 0.5),
                new ShortReadSite("G1|PAS3", 300, 3, 1.0)
        );
        List<CoverageSegment> segments = Arrays.asList(
                new CoverageSegment(0, 90, 150, "S1"),
                new CoverageSegment(1, 151, 250, "S2"),
                new CoverageSegment(2, 251, 350, "S3")
        );
        return new ShortReadGroup(
                "G1|TE001",
                "G1",
                "GENE1",
                "protein_coding",
                "chr1",
                Strand.PLUS,
                "TANDEM_3UTR",
                90,
                sites,
                segments
        );
    }


    private static ShortReadGroup fourSiteGroup() {
        List<ShortReadSite> sites = Arrays.asList(
                new ShortReadSite("G2|PAS1", 100, 1, 0.0),
                new ShortReadSite("G2|PAS2", 200, 2, 1.0 / 3.0),
                new ShortReadSite("G2|PAS3", 300, 3, 2.0 / 3.0),
                new ShortReadSite("G2|PAS4", 400, 4, 1.0)
        );
        List<CoverageSegment> segments = Arrays.asList(
                new CoverageSegment(0, 90, 150, "S1"),
                new CoverageSegment(1, 151, 250, "S2"),
                new CoverageSegment(2, 251, 350, "S3"),
                new CoverageSegment(3, 351, 450, "S4")
        );
        return new ShortReadGroup(
                "G2|TE001",
                "G2",
                "GENE2",
                "protein_coding",
                "chr1",
                Strand.PLUS,
                "TANDEM_3UTR",
                90,
                sites,
                segments
        );
    }

    private static ShortReadQuantification quant4(
            ShortReadGroup group,
            int pas1,
            int pas2,
            int pas3,
            int pas4
    ) {
        int total = pas1 + pas2 + pas3 + pas4;
        int[] counts = new int[]{pas1, pas2, pas3, pas4};
        double[] pau = new double[]{
                pas1 / (double) total,
                pas2 / (double) total,
                pas3 / (double) total,
                pas4 / (double) total
        };
        return new ShortReadQuantification(
                group,
                new long[]{0, 0, 0, 0},
                new double[]{0, 0, 0, 0},
                new double[]{0, 0, 0, 0},
                new double[]{pas1, pas2, pas3, pas4},
                pau,
                counts,
                total
        );
    }

    private static ShortReadQuantification quant(
            ShortReadGroup group,
            int pas1,
            int pas2,
            int pas3,
            int total
    ) {
        int[] counts = new int[]{pas1, pas2, pas3};
        double[] pau = new double[3];
        if (total > 0) {
            pau[0] = pas1 / (double) total;
            pau[1] = pas2 / (double) total;
            pau[2] = pas3 / (double) total;
        }
        return new ShortReadQuantification(
                group,
                new long[]{0, 0, 0},
                new double[]{0, 0, 0},
                new double[]{0, 0, 0},
                new double[]{pas1, pas2, pas3},
                pau,
                counts,
                total
        );
    }

    private static List<ShortReadQuantification> singleton(
            ShortReadQuantification quantification
    ) {
        return Arrays.asList(quantification);
    }
}
