package org.apaura.stats;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;

class IsotonicRegressionTest {

    @Test
    void decreasingRegressionPoolsAdjacentViolation() {
        double[] values = {3.0, 1.0, 2.0};
        double[] weights = {1.0, 1.0, 1.0};
        double[] expected = {3.0, 1.5, 1.5};
        assertArrayEquals(expected, IsotonicRegression.decreasing(values, weights), 1.0E-12);
    }
}
