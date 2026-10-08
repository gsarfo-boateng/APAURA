package org.apaura.analysis;

import org.apaura.model.CoverageSegment;
import org.apaura.model.Sample;
import org.apaura.model.ShortReadGroup;
import org.apaura.model.ShortReadQuantification;
import org.apaura.model.ShortReadSite;
import org.apaura.model.Strand;
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
