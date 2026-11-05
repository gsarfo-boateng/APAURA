package org.apaura.analysis;

/**
 * Shared mathematics for APA usage distributions.
 *
 * <p>Usage-change magnitude and positional direction are intentionally kept
 * separate. Total-variation distance (TVD) and Jensen-Shannon divergence (JSD)
 * summarize how much a PAS-usage distribution changed, whereas weighted
 * position index (WPI) summarizes where that change moved along transcript
 * orientation.</p>
 */
public final class ApaUsageMetrics {

    private ApaUsageMetrics() {
        // Utility class.
    }

    public static double weightedPositionIndex(
            double[] usage,
            double[] normalizedPositions
    ) {
        requireSameLength(usage, normalizedPositions);
        if (usage.length == 0) {
            return Double.NaN;
        }

        double total = positiveSum(usage);
        if (total <= 0.0) {
            return Double.NaN;
        }

        double weighted = 0.0;
        for (int index = 0; index < usage.length; index++) {
            double value = positiveFinite(usage[index]);
            weighted += (value / total) * normalizedPositions[index];
        }
        return weighted;
    }

    public static double totalVariationDistance(double[] left, double[] right) {
        requireSameLength(left, right);
        if (left.length == 0) {
            return Double.NaN;
        }

        double[] p = normalizedCopy(left);
        double[] q = normalizedCopy(right);
        if (p == null || q == null) {
            return Double.NaN;
        }

        double distance = 0.0;
        for (int index = 0; index < p.length; index++) {
            distance += Math.abs(p[index] - q[index]);
        }
        return 0.5 * distance;
    }

    /** Jensen-Shannon divergence in bits, bounded to [0, 1]. */
    public static double jensenShannonDivergence(double[] left, double[] right) {
        requireSameLength(left, right);
        if (left.length == 0) {
            return Double.NaN;
        }

        double[] p = normalizedCopy(left);
        double[] q = normalizedCopy(right);
        if (p == null || q == null) {
            return Double.NaN;
        }

        double divergence = 0.0;
        for (int index = 0; index < p.length; index++) {
            double midpoint = 0.5 * (p[index] + q[index]);
            if (p[index] > 0.0) {
                divergence += 0.5 * p[index] * log2(p[index] / midpoint);
            }
            if (q[index] > 0.0) {
                divergence += 0.5 * q[index] * log2(q[index] / midpoint);
            }
        }

        return Math.max(0.0, Math.min(1.0, divergence));
    }

    public static String observedDirection(double deltaWpi, double threshold) {
        if (!Double.isFinite(deltaWpi)) {
            return "not_applicable";
        }
        if (deltaWpi >= threshold) {
            return "lengthening";
        }
        if (deltaWpi <= -threshold) {
            return "shortening";
        }
        return "redistribution";
    }

    public static String supportedDirection(
            boolean observedChange,
            double deltaWpi,
            double threshold
    ) {
        if (!observedChange) {
            return "stable";
        }
        return observedDirection(deltaWpi, threshold);
    }

    static void addInPlace(double[] target, double[] values) {
        requireSameLength(target, values);
        for (int index = 0; index < target.length; index++) {
            target[index] += values[index];
        }
    }

    static void averageInPlace(double[] values, int count) {
        if (count <= 0) {
            return;
        }
        for (int index = 0; index < values.length; index++) {
            values[index] /= count;
        }
    }

    static double mean(java.util.List<Double> values) {
        if (values.isEmpty()) {
            return Double.NaN;
        }
        double sum = 0.0;
        for (double value : values) {
            sum += value;
        }
        return sum / values.size();
    }

    static double sampleSd(java.util.List<Double> values) {
        if (values.size() < 2) {
            return Double.NaN;
        }
        double mean = mean(values);
        double ss = 0.0;
        for (double value : values) {
            double delta = value - mean;
            ss += delta * delta;
        }
        return Math.sqrt(ss / (values.size() - 1.0));
    }

    private static double[] normalizedCopy(double[] values) {
        double total = positiveSum(values);
        if (total <= 0.0) {
            return null;
        }

        double[] normalized = new double[values.length];
        for (int index = 0; index < values.length; index++) {
            normalized[index] = positiveFinite(values[index]) / total;
        }
        return normalized;
    }

    private static double positiveSum(double[] values) {
        double total = 0.0;
        for (double value : values) {
            total += positiveFinite(value);
        }
        return total;
    }

    private static double positiveFinite(double value) {
        return Double.isFinite(value) && value > 0.0 ? value : 0.0;
    }

    private static void requireSameLength(double[] left, double[] right) {
        if (left.length != right.length) {
            throw new IllegalArgumentException(
                    "APA usage vectors must have the same length"
            );
        }
    }

    private static double log2(double value) {
        return Math.log(value) / Math.log(2.0);
    }
}
