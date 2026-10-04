package org.apaura.stats;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;

class MultipleTestingTest {

    @Test
    void benjaminiHochbergReturnsExpectedAdjustedValues() {
        double[] p = {0.01, 0.04, 0.03};
        double[] expected = {0.03, 0.04, 0.04};
        assertArrayEquals(expected, MultipleTesting.bh(p), 1.0E-12);
    }
}
