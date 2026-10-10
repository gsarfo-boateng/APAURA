package org.apaura.stats;

import org.junit.jupiter.api.Test;

import java.util.LinkedHashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

class TieAdjustedKsTest {

    @Test
    void separatedPasDistributionsHaveMaximumDistanceOne() {
        Map<Integer, Integer> control = new LinkedHashMap<Integer, Integer>();
        Map<Integer, Integer> treated = new LinkedHashMap<Integer, Integer>();
        control.put(100, 20);
        treated.put(200, 20);

        TieAdjustedKs.Result result = TieAdjustedKs.test(control, treated);

        assertEquals(1.0, result.d, 1.0E-12);
        assertEquals(1.0, result.signedD, 1.0E-12);
        assertEquals(100, result.maxDiffAt0);
        assertTrue(result.pvalue < 0.05);
    }
}
