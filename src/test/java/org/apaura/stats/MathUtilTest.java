package org.apaura.stats;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class MathUtilTest {

    @Test
    void logGammaMatchesFactorialIdentity() {
        assertEquals(Math.log(24.0), MathUtil.logGamma(5.0), 1.0E-10);
    }

    @Test
    void logChooseMatchesKnownCombination() {
        assertEquals(Math.log(10.0), MathUtil.logChoose(5, 2), 1.0E-10);
    }
}
