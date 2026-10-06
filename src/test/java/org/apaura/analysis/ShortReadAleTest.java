package org.apaura.analysis;

import org.apaura.model.CoverageSegment;
import org.apaura.model.Interval;
import org.apaura.model.Sample;
import org.apaura.model.Strand;
import org.apaura.model.TranscriptModel;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

class ShortReadAleTest {

    @Test
    void buildsAllAlternativeTerminalBranchesForMixedGene() {
        List<TranscriptModel> transcripts = Arrays.asList(
                transcript("T1", 0, 100, 100, 200),
                transcript("T2", 0, 100, 100, 300),
                transcript("T3", 0, 400, 400, 500),
                transcript("T4", 0, 600, 600, 700)
        );

        List<ArchitectureCatalog.Entry> architecture =
                ArchitectureCatalog.build(transcripts, 5, 0);
        assertEquals(1, architecture.size());
        assertEquals("mixed", architecture.get(0).architectureClass);
        assertEquals("TUTR+ALE", architecture.get(0).components);
        assertEquals(3, architecture.get(0).terminalBranches);
        assertEquals(1, architecture.get(0).tutrBranches);

        List<ShortReadAle.Group> groups = ShortReadAle.build(transcripts, 5, 0, 10);
        assertEquals(1, groups.size());
        assertEquals(3, groups.get(0).branches.size());
        assertEquals(2, groups.get(0).branches.get(0).annotatedPas0.size());
        assertEquals(1, groups.get(0).branches.get(1).annotatedPas0.size());
        assertEquals(1, groups.get(0).branches.get(2).annotatedPas0.size());

        // ALE evidence for a TUTR-capable branch must stop at its proximal PAS.
        // Sequence unique only to the distal TUTR isoform must not inflate ALE usage.
        assertEquals(100, groups.get(0).branches.get(0).availableUniqueBases);
        assertEquals(100, groups.get(0).branches.get(0).evidenceBases);
        assertEquals(1, groups.get(0).branches.get(0).uniqueRegions.size());
        assertEquals(100, groups.get(0).branches.get(0).uniqueRegions.get(0).start0);
        assertEquals(200, groups.get(0).branches.get(0).uniqueRegions.get(0).end0);
    }

    @Test
    void detectsReplicatedAleSwitchWithBetaBinomial() {
        ShortReadAle.Group group = twoBranchGroup();
        List<Sample> samples = new ArrayList<Sample>();
        Map<String, List<ShortReadAle.Quantification>> bySample =
                new LinkedHashMap<String, List<ShortReadAle.Quantification>>();

        for (int i = 1; i <= 3; i++) {
            String name = "WT" + i;
            samples.add(new Sample(name, "control", name + ".bam"));
            bySample.put(name, singleton(new ShortReadAle.Quantification(group, new int[]{90, 10})));
        }
        for (int i = 1; i <= 3; i++) {
            String name = "MUT" + i;
            samples.add(new Sample(name, "case", name + ".bam"));
            bySample.put(name, singleton(new ShortReadAle.Quantification(group, new int[]{10, 90})));
        }

        ShortReadAle.Results results = ShortReadAle.differential(
                samples,
                bySample,
                "control",
                "case",
                10,
                3,
                0.10,
                0.10,
                0.05,
                0.10
        );

        assertEquals(2, results.branches.size());
        assertEquals(1, results.genes.size());
        assertEquals("ale-beta-binomial-LRT", results.branches.get(0).method);
        assertTrue(Double.isFinite(results.branches.get(0).q));
        assertTrue(results.genes.get(0).aleChange);
        assertEquals("distal_switch", results.genes.get(0).aleDirection);
        assertTrue(results.genes.get(0).switchIndex > 0.5);
    }

    @Test
    void unchangedBranchUsageIsStableEvenWhenGeneHasMixedArchitecture() {
        ShortReadAle.Group group = twoBranchGroup();
        List<Sample> samples = new ArrayList<Sample>();
        Map<String, List<ShortReadAle.Quantification>> bySample =
                new LinkedHashMap<String, List<ShortReadAle.Quantification>>();

        for (int i = 1; i <= 3; i++) {
            String name = "WT" + i;
            samples.add(new Sample(name, "control", name + ".bam"));
            bySample.put(name, singleton(new ShortReadAle.Quantification(group, new int[]{80, 20})));
        }
        for (int i = 1; i <= 3; i++) {
            String name = "MUT" + i;
            samples.add(new Sample(name, "case", name + ".bam"));
            bySample.put(name, singleton(new ShortReadAle.Quantification(group, new int[]{80, 20})));
        }

        ShortReadAle.Results results = ShortReadAle.differential(
                samples, bySample, "control", "case",
                10, 3, 0.10, 0.10, 0.05, 0.10
        );

        assertFalse(results.genes.get(0).aleChange);
        assertEquals(0.0, results.genes.get(0).switchIndex, 1.0E-12);
        assertEquals("stable", results.genes.get(0).observedDirection);
        assertEquals("stable", results.genes.get(0).aleDirection);
    }

    @Test
    void sparseReplicatedAleIsNotRescuedByPooledFisher() {
        ShortReadAle.Group group = twoBranchGroup();
        List<Sample> samples = new ArrayList<Sample>();
        Map<String, List<ShortReadAle.Quantification>> bySample =
                new LinkedHashMap<String, List<ShortReadAle.Quantification>>();

        for (int i = 1; i <= 10; i++) {
            String name = "WT" + i;
            samples.add(new Sample(name, "control", name + ".bam"));
            bySample.put(name, singleton(new ShortReadAle.Quantification(group, new int[]{0, 0})));
        }
        for (int i = 1; i <= 10; i++) {
            String name = "MUT" + i;
            samples.add(new Sample(name, "case", name + ".bam"));
            bySample.put(name, singleton(new ShortReadAle.Quantification(group, new int[]{0, 0})));
        }
        bySample.put("WT1", singleton(new ShortReadAle.Quantification(group, new int[]{90, 10})));
        bySample.put("MUT1", singleton(new ShortReadAle.Quantification(group, new int[]{10, 90})));
        bySample.put("MUT2", singleton(new ShortReadAle.Quantification(group, new int[]{10, 90})));

        ShortReadAle.Results results = ShortReadAle.differential(
                samples,
                bySample,
                "control",
                "case",
                10,
                3,
                0.10,
                0.10,
                0.05,
                0.10
        );

        assertEquals("insufficient_replicate_support", results.genes.get(0).status);
        assertFalse(results.genes.get(0).aleChange);
        assertTrue(Double.isNaN(results.branches.get(0).p));
        assertTrue(Double.isNaN(results.branches.get(0).q));
    }

    private static TranscriptModel transcript(
            String id,
            int upstreamStart,
            int upstreamEnd,
            int terminalStart,
            int terminalEnd
    ) {
        TranscriptModel tx = new TranscriptModel(
                id,
                "G1",
                "GENE1",
                "protein_coding",
                "chr1",
                Strand.PLUS
        );
        tx.addExon(upstreamStart, upstreamEnd);
        tx.addExon(terminalStart, terminalEnd);
        tx.addCds(upstreamStart, Math.min(upstreamEnd, upstreamStart + 50));
        return tx;
    }

    private static ShortReadAle.Group twoBranchGroup() {
        ShortReadAle.Branch b1 = new ShortReadAle.Branch(
                "G1|ALEB001", "G1", "GENE1", "protein_coding", "chr1",
                Strand.PLUS, 1, 100, 200, Arrays.asList(200),
                Arrays.asList(new CoverageSegment(0, 100, 200, "ALE_UNIQUE")),
                100, 100, 0.0
        );
        ShortReadAle.Branch b2 = new ShortReadAle.Branch(
                "G1|ALEB002", "G1", "GENE1", "protein_coding", "chr1",
                Strand.PLUS, 2, 300, 400, Arrays.asList(400),
                Arrays.asList(new CoverageSegment(0, 300, 400, "ALE_UNIQUE")),
                100, 100, 1.0
        );
        return new ShortReadAle.Group(
                "G1|ALE", "G1", "GENE1", "protein_coding", "chr1",
                Strand.PLUS, Arrays.asList(b1, b2)
        );
    }

    private static List<ShortReadAle.Quantification> singleton(
            ShortReadAle.Quantification q
    ) {
        return Arrays.asList(q);
    }
}
