package org.apaura.stats;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

class DirichletMultinomialTest {

    @Test
    void fixedPrecisionLrtIsNearNullForIdenticalComposition() {
        int[][] counts = {
                {80, 10, 10},
                {78, 12, 10},
                {82, 8, 10},
                {80, 10, 10},
                {78, 12, 10},
                {82, 8, 10}
        };
        int[] groups = {0, 0, 0, 1, 1, 1};

        DirichletMultinomial.TestResult result =
                DirichletMultinomial.testTwoGroup(counts, groups, 25.0);

        assertTrue(result.converged);
        assertEquals(2, result.degreesOfFreedom);
        assertTrue(result.likelihoodRatio >= 0.0);
        assertTrue(result.likelihoodRatio < 1.0E-8);
        assertTrue(result.pvalue > 0.999);
    }

    @Test
    void fixedPrecisionLrtDetectsStrongCompositionShift() {
        int[][] counts = {
                {80, 10, 10},
                {78, 12, 10},
                {82, 8, 10},
                {10, 10, 80},
                {12, 8, 80},
                {8, 12, 80}
        };
        int[] groups = {0, 0, 0, 1, 1, 1};

        DirichletMultinomial.TestResult result =
                DirichletMultinomial.testTwoGroup(counts, groups, 25.0);

        assertTrue(result.converged);
        assertEquals(2, result.degreesOfFreedom);
        assertTrue(result.fullLogLikelihood > result.nullLogLikelihood);
        assertTrue(result.likelihoodRatio > 10.0);
        assertTrue(result.pvalue < 0.01);
    }

    @Test
    void logLikelihoodIsFiniteForValidComposition() {
        int[][] counts = {
                {30, 20, 50},
                {25, 25, 50},
                {35, 15, 50}
        };

        double ll = DirichletMultinomial.logLikelihood(
                counts,
                new double[]{0.30, 0.20, 0.50},
                20.0
        );

        assertTrue(Double.isFinite(ll));
    }
}
