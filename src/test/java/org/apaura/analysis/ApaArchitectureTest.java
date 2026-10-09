package org.apaura.analysis;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class ApaArchitectureTest {

    @Test
    void classifiesCurrentTerminalArchitecturesConsistently() {
        assertEquals(
                "none",
                ApaArchitecture.classify(false, false)
        );

        assertEquals(
                "TUTR",
                ApaArchitecture.classify(true, false)
        );

        assertEquals(
                "ALE",
                ApaArchitecture.classify(false, true)
        );

        assertEquals(
                "mixed",
                ApaArchitecture.classify(true, true)
        );
    }

    @Test
    void reportsCurrentArchitectureComponentsConsistently() {
        assertEquals(
                "none",
                ApaArchitecture.components(false, false)
        );

        assertEquals(
                "TUTR",
                ApaArchitecture.components(true, false)
        );

        assertEquals(
                "ALE",
                ApaArchitecture.components(false, true)
        );

        assertEquals(
                "TUTR+ALE",
                ApaArchitecture.components(true, true)
        );
    }
}