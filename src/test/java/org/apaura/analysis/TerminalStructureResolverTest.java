package org.apaura.analysis;

import org.apaura.model.Sample;
import org.apaura.model.SiteCluster;
import org.apaura.model.Strand;
import org.apaura.model.TranscriptModel;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;

class TerminalStructureResolverTest {

    @Test
    void sharedPhysicalPasConnectsRawAnchorsIntoOneTerminalStructure() {
        TranscriptModel first = transcript(
                "T1",
                100,
                500
        );

        TranscriptModel second = transcript(
                "T2",
                140,
                500
        );

        TranscriptModel third = transcript(
                "T3",
                100,
                600
        );

        SiteCluster proximal = site(
                "PAS1",
                499,
                140
        );

        SiteCluster distal = site(
                "PAS2",
                599,
                100
        );

        TerminalStructureResolver.assign(
                Arrays.asList(
                        first,
                        second,
                        third
                ),
                Arrays.asList(
                        proximal,
                        distal
                ),
                0,
                24
        );

        assertEquals(
                proximal.terminalStructureAnchor0,
                distal.terminalStructureAnchor0
        );

        assertEquals(
                140,
                proximal.terminalStructureAnchor0
        );

        ArchitectureAnalysis.annotateSiteRoles(
                Arrays.asList(
                        proximal,
                        distal
                )
        );

        assertEquals(1, proximal.terminalBranchRank);
        assertEquals(1, distal.terminalBranchRank);

        assertEquals("TUTR", proximal.architectureRole);
        assertEquals("TUTR", distal.architectureRole);
    }


    @Test
    void retainedTerminalPasWithoutResolvedStructureForcesUnresolvedArchitecture() {
        SiteCluster resolved = site(
                "PAS1",
                499,
                100
        );
        resolved.terminalStructureAnchor0 = 100;
        resolved.region = "3UTR";

        SiteCluster unresolved = site(
                "PAS2",
                599,
                140
        );
        unresolved.terminalStructureAnchor0 = -1;
        unresolved.region = "3UTR";

        resolved.counts.put("C1", 10);
        resolved.counts.put("T1", 10);
        unresolved.counts.put("C1", 10);
        unresolved.counts.put("T1", 10);

        Map<String, ArchitectureAnalysis.Metrics> metrics =
                ArchitectureAnalysis.analyze(
                        Arrays.asList(
                                new Sample("C1", "control", "C1.bam"),
                                new Sample("T1", "case", "T1.bam")
                        ),
                        Arrays.asList(resolved, unresolved),
                        "control",
                        "case"
                );

        ArchitectureAnalysis.Metrics gene = metrics.get("G1");
        assertEquals("UNRESOLVED", gene.architectureClass);
        assertEquals(1, gene.nTerminalBranches);
        assertEquals(0, gene.nTutrBranches);
        assertEquals("terminal", resolved.architectureRole);
        assertEquals("other", unresolved.architectureRole);
    }

    private static TranscriptModel transcript(
            String transcriptId,
            int terminalAnchor0,
            int transcriptEnd0
    ) {
        TranscriptModel transcript = new TranscriptModel(
                transcriptId,
                "G1",
                "GENE1",
                "protein_coding",
                "chr1",
                Strand.PLUS
        );

        transcript.addExon(
                terminalAnchor0,
                transcriptEnd0
        );

        transcript.addCds(
                terminalAnchor0,
                Math.min(
                        terminalAnchor0 + 100,
                        transcriptEnd0
                )
        );

        return transcript;
    }

    private static SiteCluster site(
            String id,
            int position0,
            int rawAnchor0
    ) {
        SiteCluster site = new SiteCluster(
                id,
                "chr1",
                Strand.PLUS,
                position0,
                position0 + 1,
                position0
        );

        site.geneId = "G1";
        site.geneName = "GENE1";
        site.biotype = "protein_coding";
        site.terminalExon = true;
        site.terminalExonAnchor0 = rawAnchor0;

        return site;
    }
}