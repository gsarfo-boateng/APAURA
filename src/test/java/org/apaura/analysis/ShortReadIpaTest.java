package org.apaura.analysis;

import org.apaura.model.Strand;
import org.apaura.model.TranscriptModel;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

class ShortReadIpaTest {

    @Test
    void buildsPlusStrandAnnotatedIpaEvent() {
        TranscriptModel full = transcript("full", "G1", Strand.PLUS);
        full.addExon(100, 200);
        full.addExon(500, 600);

        TranscriptModel ipa = transcript("ipa", "G1", Strand.PLUS);
        ipa.addExon(100, 200);
        ipa.addExon(300, 350);

        List<ShortReadIpa.Group> groups = ShortReadIpa.build(
                Arrays.asList(full, ipa),
                24,
                100,
                30
        );

        assertEquals(1, groups.size());
        ShortReadIpa.Group group = groups.get(0);
        assertEquals(349, group.representativePas0);
        assertEquals(150, group.commonWindow.start0);
        assertEquals(200, group.commonWindow.end0);
        assertEquals(300, group.ipaSpecificWindow.start0);
        assertEquals(350, group.ipaSpecificWindow.end0);
        assertEquals(500, group.continuationWindow.start0);
        assertEquals(550, group.continuationWindow.end0);
    }

    @Test
    void buildsMinusStrandAnnotatedIpaEvent() {
        TranscriptModel full = transcript("full", "G2", Strand.MINUS);
        full.addExon(100, 200);
        full.addExon(500, 600);

        TranscriptModel ipa = transcript("ipa", "G2", Strand.MINUS);
        ipa.addExon(500, 600);
        ipa.addExon(350, 400);

        List<ShortReadIpa.Group> groups = ShortReadIpa.build(
                Arrays.asList(full, ipa),
                24,
                100,
                30
        );

        assertEquals(1, groups.size());
        ShortReadIpa.Group group = groups.get(0);
        assertEquals(350, group.representativePas0);
        assertEquals(500, group.commonWindow.start0);
        assertEquals(550, group.commonWindow.end0);
        assertEquals(350, group.ipaSpecificWindow.start0);
        assertEquals(400, group.ipaSpecificWindow.end0);
        assertEquals(150, group.continuationWindow.start0);
        assertEquals(200, group.continuationWindow.end0);
    }


    @Test
    void clustersNearbyIpaPasButKeepsDistantSitesSeparate() {
        TranscriptModel full = transcript("full", "G4", Strand.PLUS);
        full.addExon(100, 200);
        full.addExon(500, 600);

        TranscriptModel ipa1 = transcript("ipa1", "G4", Strand.PLUS);
        ipa1.addExon(100, 200);
        ipa1.addExon(300, 350);

        TranscriptModel ipa2 = transcript("ipa2", "G4", Strand.PLUS);
        ipa2.addExon(100, 200);
        ipa2.addExon(320, 360);

        TranscriptModel ipa3 = transcript("ipa3", "G4", Strand.PLUS);
        ipa3.addExon(100, 200);
        ipa3.addExon(400, 430);

        List<ShortReadIpa.Group> groups = ShortReadIpa.build(
                Arrays.asList(full, ipa1, ipa2, ipa3),
                24,
                100,
                30
        );

        assertEquals(2, groups.size());
        assertEquals(Arrays.asList(349, 359), groups.get(0).annotatedPas0);
        assertEquals(359, groups.get(0).representativePas0);
        assertEquals(Arrays.asList(429), groups.get(1).annotatedPas0);
        assertEquals(429, groups.get(1).representativePas0);
    }

    @Test
    void convertsCoverageDropToIpaUsage() {
        TranscriptModel full = transcript("full", "G3", Strand.PLUS);
        full.addExon(100, 200);
        full.addExon(500, 600);

        TranscriptModel ipa = transcript("ipa", "G3", Strand.PLUS);
        ipa.addExon(100, 200);
        ipa.addExon(300, 350);

        ShortReadIpa.Group group = ShortReadIpa.build(
                Arrays.asList(full, ipa),
                24,
                100,
                30
        ).get(0);

        ShortReadIpa.Quantification quantification =
                ShortReadIpa.quantify(group, 100, 60, 40);

        assertEquals(60, quantification.ipaSpecificReads);
        assertEquals(0.60, quantification.ipaUsage, 1.0E-12);
    }

    private static TranscriptModel transcript(
            String transcriptId,
            String geneId,
            Strand strand
    ) {
        return new TranscriptModel(
                transcriptId,
                geneId,
                geneId,
                "protein_coding",
                "chr1",
                strand
        );
    }
}
