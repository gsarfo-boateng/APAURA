package org.apaura.analysis;

import org.apaura.model.Sample;
import org.apaura.model.SiteCluster;
import org.apaura.model.Strand;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

class DifferentialApaTest {

    @Test
    void retainsVectorMagnitudeAndJsdForMultiPasRedistribution() {
        List<Sample> samples = Arrays.asList(
                new Sample("C1", "control", "C1.bam"),
                new Sample("T1", "case", "T1.bam")
        );

        List<SiteCluster> sites = fourSites();
        setCounts(sites, "C1", new int[]{250, 250, 250, 250});
        setCounts(sites, "T1", new int[]{190, 310, 190, 310});

        DifferentialApa.Results results = DifferentialApa.run(
                samples,
                sites,
                "control",
                "case",
                1,
                0.10,
                0.05,
                0.10
        );

        DifferentialApa.GeneResult gene = results.genes.get(0);
        assertEquals(0.12, gene.apd, 1.0E-12);
        assertTrue(gene.usageJsd > 0.0);
        assertEquals(0.04, gene.deltaWpi, 1.0E-12);
        assertTrue(gene.apaChange);
        assertEquals("no_directional_shift", gene.apaDirection);
    }

    @Test
    void fisherFallbackIsLimitedToGenuineOneVsOne() {
        List<Sample> samples = Arrays.asList(
                new Sample("C1", "control", "C1.bam"),
                new Sample("T1", "case", "T1.bam"),
                new Sample("T2", "case", "T2.bam")
        );

        List<SiteCluster> sites = twoSites();
        setCounts(sites, "C1", new int[]{90, 10});
        setCounts(sites, "T1", new int[]{10, 90});
        setCounts(sites, "T2", new int[]{12, 88});

        DifferentialApa.Results results = DifferentialApa.run(
                samples,
                sites,
                "control",
                "case",
                1,
                0.10,
                0.05,
                0.10
        );

        for (DifferentialApa.SiteResult site : results.sites) {
            assertEquals(1, site.controlReplicates);
            assertEquals(2, site.caseReplicates);
            assertEquals("not_tested", site.method);
            assertEquals("insufficient_replicate_support", site.status);
            assertTrue(Double.isNaN(site.p));
            assertTrue(Double.isNaN(site.q));
        }
    }

    @Test
    void genuineOneVsOneUsesFisherExact() {
        List<Sample> samples = Arrays.asList(
                new Sample("C1", "control", "C1.bam"),
                new Sample("T1", "case", "T1.bam")
        );

        List<SiteCluster> sites = twoSites();
        setCounts(sites, "C1", new int[]{90, 10});
        setCounts(sites, "T1", new int[]{10, 90});

        DifferentialApa.Results results = DifferentialApa.run(
                samples,
                sites,
                "control",
                "case",
                1,
                0.10,
                0.05,
                0.10
        );

        for (DifferentialApa.SiteResult site : results.sites) {
            assertEquals("Fisher-exact-unreplicated", site.method);
            assertEquals("ok", site.status);
            assertTrue(Double.isFinite(site.p));
            assertTrue(Double.isFinite(site.q));
        }
    }


    @Test
    void changedMechanismRequiresSignificantSiteWithinThatMechanism() {
        List<Sample> samples = Arrays.asList(
                new Sample("C1", "control", "C1.bam"),
                new Sample("T1", "case", "T1.bam")
        );

        SiteCluster tutrProximal = site("PAS1", 100, 1, 0.0);
        SiteCluster tutrDistal = site("PAS2", 200, 2, 1.0 / 3.0);

        SiteCluster internal = new SiteCluster(
                "PAS3", "chr1", Strand.PLUS, 300, 300, 300
        );
        internal.geneId = "G1";
        internal.geneName = "GENE1";
        internal.biotype = "protein_coding";
        internal.region = "EXON_CODING";
        internal.siteRank = 3;
        internal.normalizedPosition = 2.0 / 3.0;

        SiteCluster downstream = new SiteCluster(
                "PAS4", "chr1", Strand.PLUS, 400, 400, 400
        );
        downstream.geneId = "G1";
        downstream.geneName = "GENE1";
        downstream.biotype = "protein_coding";
        downstream.region = "DOWNSTREAM_3P";
        downstream.siteRank = 4;
        downstream.normalizedPosition = 1.0;

        List<SiteCluster> sites = Arrays.asList(
                tutrProximal,
                tutrDistal,
                internal,
                downstream
        );

        setCounts(sites, "C1", new int[]{6, 4, 100, 0});
        setCounts(sites, "T1", new int[]{4, 6, 0, 100});

        ArchitectureAnalysis.annotateSiteRoles(sites);

        DifferentialApa.Results differential = DifferentialApa.run(
                samples,
                sites,
                "control",
                "case",
                1,
                0.10,
                0.05,
                0.10
        );

        ArchitectureAnalysis.enrichGeneResults(
                differential.genes,
                ArchitectureAnalysis.analyze(
                        samples,
                        sites,
                        "control",
                        "case"
                ),
                0.05,
                0.10,
                0.10,
                0.10,
                0.10
        );

        DifferentialApa.GeneResult gene = differential.genes.get(0);

        assertTrue(gene.tutrShiftMagnitude >= 0.10);
        assertTrue(gene.tutrMinFdr > 0.05);
        assertTrue(gene.internalExonMinFdr <= 0.05);
        assertEquals("stable", gene.tutrDirection);
        assertEquals("internal_exon", gene.changedMechanisms);
        assertEquals("internal_exon:loss", gene.apaEvent);
        assertTrue(gene.apaChange);
    }

    @Test
    void unresolvedArchitectureCanRetainSupportedLocalTutrSignal() {
        List<Sample> samples = Arrays.asList(
                new Sample("C1", "control", "C1.bam"),
                new Sample("T1", "case", "T1.bam")
        );

        SiteCluster proximal = site("PAS1", 100, 1, 0.0);
        SiteCluster distal = site("PAS2", 200, 2, 0.5);

        SiteCluster unresolved = site("PAS3", 300, 3, 1.0);
        unresolved.terminalStructureAnchor0 = -1;
        unresolved.region = "3UTR";

        List<SiteCluster> sites = Arrays.asList(
                proximal, distal, unresolved
        );

        setCounts(sites, "C1", new int[]{90, 10, 10});
        setCounts(sites, "T1", new int[]{10, 90, 10});

        ArchitectureAnalysis.annotateSiteRoles(sites);

        DifferentialApa.Results differential = DifferentialApa.run(
                samples,
                sites,
                "control",
                "case",
                1,
                0.10,
                0.05,
                0.10
        );

        ArchitectureAnalysis.enrichGeneResults(
                differential.genes,
                ArchitectureAnalysis.analyze(
                        samples,
                        sites,
                        "control",
                        "case"
                ),
                0.05,
                0.10,
                0.10,
                0.10,
                0.10
        );

        DifferentialApa.GeneResult gene = differential.genes.get(0);

        assertEquals("UNRESOLVED", gene.apaType);
        assertEquals("TUTR", gene.architectureComponents);
        assertTrue(gene.tutrMinFdr <= 0.05);
        assertEquals("TUTR", gene.changedMechanisms);
        assertEquals("lengthening", gene.tutrDirection);
        assertTrue(gene.apaChange);
    }

    private static List<SiteCluster> fourSites() {
        List<SiteCluster> sites = new ArrayList<SiteCluster>();
        sites.add(site("PAS1", 100, 1, 0.0));
        sites.add(site("PAS2", 200, 2, 1.0 / 3.0));
        sites.add(site("PAS3", 300, 3, 2.0 / 3.0));
        sites.add(site("PAS4", 400, 4, 1.0));
        return sites;
    }

    private static List<SiteCluster> twoSites() {
        return Arrays.asList(
                site("PAS1", 100, 1, 0.0),
                site("PAS2", 200, 2, 1.0)
        );
    }

    private static SiteCluster site(
            String id,
            int position0,
            int rank,
            double normalizedPosition
    ) {
        SiteCluster site = new SiteCluster(
                id,
                "chr1",
                Strand.PLUS,
                position0,
                position0,
                position0
        );
        site.geneId = "G1";
        site.geneName = "GENE1";
        site.biotype = "protein_coding";
        site.siteRank = rank;
        site.normalizedPosition = normalizedPosition;
        site.terminalExon = true;
        site.terminalExonAnchor0 = 50;
        site.terminalStructureAnchor0 = 50;
        return site;
    }

    private static void setCounts(
            List<SiteCluster> sites,
            String sample,
            int[] counts
    ) {
        for (int index = 0; index < sites.size(); index++) {
            sites.get(index).counts.put(sample, counts[index]);
        }
    }
}
