package org.apaura.analysis;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

class ApaArchitectureTest {

    @Test
    void classifiesTerminalArchitectureFromStructureAndPasCounts() {
        assertEquals(
                "UNRESOLVED",
                ApaArchitecture.classify(0, 0, 0)
        );

        assertEquals(
                "SINGLE_PAS",
                ApaArchitecture.classify(1, 0, 1)
        );

        assertEquals(
                "TUTR",
                ApaArchitecture.classify(1, 1, 2)
        );

        assertEquals(
                "ALE",
                ApaArchitecture.classify(2, 0, 2)
        );

        assertEquals(
                "mixed",
                ApaArchitecture.classify(2, 1, 3)
        );
    }

    @Test
    void reportsOnlyAlternativeArchitectureComponents() {
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

    @Test
    void distinguishesAlternativeArchitectureFromSinglePasStructure() {
        assertFalse(
                ApaArchitecture.isAlternativeTerminalArchitecture("SINGLE_PAS")
        );

        assertFalse(
                ApaArchitecture.isAlternativeTerminalArchitecture("UNRESOLVED")
        );

        assertTrue(
                ApaArchitecture.isAlternativeTerminalArchitecture("TUTR")
        );

        assertTrue(
                ApaArchitecture.isAlternativeTerminalArchitecture("ALE")
        );

        assertTrue(
                ApaArchitecture.isAlternativeTerminalArchitecture("mixed")
        );
    }
}