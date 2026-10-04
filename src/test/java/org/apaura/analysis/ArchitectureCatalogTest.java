package org.apaura.analysis;

import org.apaura.model.Strand;
import org.apaura.model.TranscriptModel;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

class ArchitectureCatalogTest {

    @Test
    void classicalIpaDoesNotCreateAleArchitecture() {
        TranscriptModel full = transcript("full", "G1", Strand.PLUS);
        full.addExon(100, 200);
        full.addExon(500, 600);

        TranscriptModel ipa = transcript("ipa", "G1", Strand.PLUS);
        ipa.addExon(100, 350);

        List<ArchitectureCatalog.Entry> entries = ArchitectureCatalog.build(
                Arrays.asList(full, ipa), 5, 24
        );

        assertEquals(1, entries.size());
        ArchitectureCatalog.Entry entry = entries.get(0);
        assertEquals("none", entry.architectureClass);
        assertEquals("none", entry.components);
        assertEquals(2, entry.transcripts);
        assertEquals(1, entry.terminalBranches);
        assertEquals(2, entry.annotatedPas);
        assertEquals(1, entry.terminalPas);
        assertEquals(0, entry.tutrBranches);
        assertEquals(1, entry.ipaCompatibleEnds);

        assertEquals(0, ShortReadAle.build(
                Arrays.asList(full, ipa), 5, 24, 10
        ).size());
        assertEquals(0, ShortReadGroupBuilder.build(
                Arrays.asList(full, ipa), 5, 24, 10
        ).size());
    }

    @Test
    void splicedAlternativeLastExonRemainsAleNotIpa() {
        TranscriptModel full = transcript("full", "G2", Strand.PLUS);
        full.addExon(100, 200);
        full.addExon(500, 600);

        TranscriptModel ale = transcript("ale", "G2", Strand.PLUS);
        ale.addExon(100, 200);
        ale.addExon(300, 350);

        List<ArchitectureCatalog.Entry> entries = ArchitectureCatalog.build(
                Arrays.asList(full, ale), 5, 24
        );

        assertEquals(1, entries.size());
        ArchitectureCatalog.Entry entry = entries.get(0);
        assertEquals("ALE", entry.architectureClass);
        assertEquals("ALE", entry.components);
        assertEquals(2, entry.terminalBranches);
        assertEquals(2, entry.annotatedPas);
        assertEquals(2, entry.terminalPas);
        assertEquals(0, entry.tutrBranches);
        assertEquals(0, entry.ipaCompatibleEnds);
        assertEquals(1, ShortReadAle.build(
                Arrays.asList(full, ale), 5, 24, 10
        ).size());
        assertEquals(0, ShortReadIpa.build(
                Arrays.asList(full, ale), 24, 100, 30
        ).size());
    }

    @Test
    void mixedTerminalArchitectureStaysMixed() {
        TranscriptModel tutr1 = transcript("tutr1", "G3", Strand.PLUS);
        tutr1.addExon(0, 100);
        tutr1.addExon(100, 200);

        TranscriptModel tutr2 = transcript("tutr2", "G3", Strand.PLUS);
        tutr2.addExon(0, 100);
        tutr2.addExon(100, 300);

        TranscriptModel ale = transcript("ale", "G3", Strand.PLUS);
        ale.addExon(0, 100);
        ale.addExon(400, 500);

        List<ArchitectureCatalog.Entry> entries = ArchitectureCatalog.build(
                Arrays.asList(tutr1, tutr2, ale), 5, 0
        );

        assertEquals(1, entries.size());
        ArchitectureCatalog.Entry entry = entries.get(0);
        assertEquals("mixed", entry.architectureClass);
        assertEquals("TUTR+ALE", entry.components);
        assertEquals(2, entry.terminalBranches);
        assertEquals(3, entry.annotatedPas);
        assertEquals(3, entry.terminalPas);
        assertEquals(1, entry.tutrBranches);
        assertEquals(0, entry.ipaCompatibleEnds);
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
