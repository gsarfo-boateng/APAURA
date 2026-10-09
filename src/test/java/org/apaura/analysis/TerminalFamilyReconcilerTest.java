package org.apaura.analysis;

import org.apaura.model.Strand;
import org.apaura.model.TranscriptModel;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertEquals;

class TerminalFamilyReconcilerTest {

    @Test
    void sharedPhysicalPasMergesDifferentRawTerminalAnchors() {
        TranscriptModel first = new TranscriptModel(
                "T1",
                "G1",
                "GENE1",
                "protein_coding",
                "chr1",
                Strand.PLUS
        );
        first.addExon(100, 500);
        first.addCds(120, 300);

        TranscriptModel second = new TranscriptModel(
                "T2",
                "G1",
                "GENE1",
                "protein_coding",
                "chr1",
                Strand.PLUS
        );
        second.addExon(140, 500);
        second.addCds(160, 300);

        TerminalFamilyReconciler.Result result =
                TerminalFamilyReconciler.reconcile(
                        Arrays.asList(first, second),
                        0,
                        24
                );

        assertEquals(2, result.initialAnchorBuckets);
        assertEquals(1, result.mergedAnchorBuckets);
        assertEquals(1, result.families.size());

        TerminalFamilyReconciler.Family family =
                result.families.get(0);

        assertEquals(2, family.sourceAnchorBuckets);
        assertEquals(2, family.transcriptCount);

        assertEquals(
                Arrays.asList(100, 140),
                family.sourceAnchors
        );

        assertEquals(
                Arrays.asList(499),
                family.pasPositions
        );
    }
}