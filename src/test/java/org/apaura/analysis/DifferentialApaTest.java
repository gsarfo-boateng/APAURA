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
