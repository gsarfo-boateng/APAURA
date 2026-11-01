package org.apaura.analysis;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class ProtocolRegistryTest {

    @Test
    public void genericUsesTranscriptThreePrimeBoundary() {
        assertEquals(
                ProtocolSpec.Boundary.TRANSCRIPT_3P,
                ProtocolRegistry.get("generic").boundary
        );
    }

    @Test
    public void quantSeqRevUsesFivePrimeBoundaryAndFlip() {
        ProtocolSpec spec = ProtocolRegistry.get("quantseq-rev");
        assertEquals(ProtocolSpec.Boundary.ALIGNMENT_5P, spec.boundary);
        assertTrue(spec.flipAlignmentStrand);
    }

    @Test
    public void aliasesResolve() {
        assertEquals("3p-seq", ProtocolRegistry.get("3pseq").name);
        assertEquals("a-seq", ProtocolRegistry.get("aseq").name);
    }

    @Test
    public void unknownProtocolFails() {
        assertThrows(
                IllegalArgumentException.class,
                () -> ProtocolRegistry.get("not-a-protocol")
        );
    }
}
