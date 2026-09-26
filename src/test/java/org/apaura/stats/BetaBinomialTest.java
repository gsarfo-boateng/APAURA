package org.apaura.stats;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

class BetaBinomialTest {

    @Test
    void detectsLargeReplicateConsistentUsageDifference() {
        int[] successes = {10, 12, 80, 82};
        int[] totals = {100, 100, 100, 100};
        int[] groups = {0, 0, 1, 1};

        BetaBinomial.TestResult result = BetaBinomial.test(successes, totals, groups);

        assertTrue(result.pvalue < 0.05);
        assertTrue(result.p0 < result.p1);
        assertTrue(result.rho > 0.0);
    }
}
