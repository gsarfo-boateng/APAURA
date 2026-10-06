package org.apaura.stats;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

class FisherExactTest {

    @Test
    void twoSidedTestIsSymmetricWhenRowsAreSwapped() {
        double first = FisherExact.twoSided(8, 2, 1, 9);
        double swapped = FisherExact.twoSided(1, 9, 8, 2);
        assertEquals(first, swapped, 1.0E-12);
        assertTrue(first < 0.01);
    }

    @Test
    void rightTailDetectsEnrichment() {
        assertTrue(FisherExact.rightTail(8, 2, 1, 9) < 0.01);
    }
}
