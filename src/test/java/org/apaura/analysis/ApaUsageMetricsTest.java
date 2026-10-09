package org.apaura.analysis;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

class ApaUsageMetricsTest {

    @Test
    void twoPasWpiEqualsDistalUsage() {
        double[] usage = {0.25, 0.75};
        double[] positions = {0.0, 1.0};

        assertEquals(
                0.75,
                ApaUsageMetrics.weightedPositionIndex(usage, positions),
                1.0E-12
        );
    }

    @Test
    void tvdAndJsdCaptureDistributionalChange() {
        double[] control = {1.0, 0.0};
        double[] treated = {0.0, 1.0};

        assertEquals(
                1.0,
                ApaUsageMetrics.totalVariationDistance(control, treated),
                1.0E-12
        );
        assertEquals(
                1.0,
                ApaUsageMetrics.jensenShannonDivergence(control, treated),
                1.0E-12
        );
    }

    @Test
    void multiPasRedistributionCanBeLargeWithoutLargeWpiShift() {
        double[] control = {0.25, 0.25, 0.25, 0.25};
        double[] treated = {0.19, 0.31, 0.19, 0.31};
        double[] positions = {0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0};

        double controlWpi = ApaUsageMetrics.weightedPositionIndex(
                control,
                positions
        );
        double treatedWpi = ApaUsageMetrics.weightedPositionIndex(
                treated,
                positions
        );
        double tvd = ApaUsageMetrics.totalVariationDistance(
                control,
                treated
        );
        double jsd = ApaUsageMetrics.jensenShannonDivergence(
                control,
                treated
        );

        assertEquals(0.12, tvd, 1.0E-12);
        assertEquals(0.04, treatedWpi - controlWpi, 1.0E-12);
        assertTrue(jsd > 0.0);
        assertEquals(
                "no_directional_shift",
                ApaUsageMetrics.observedDirection(
                        treatedWpi - controlWpi,
                        0.10
                )
        );
    }
}
