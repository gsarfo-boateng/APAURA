package org.apaura.stats;

import java.util.Map;
import java.util.TreeSet;

/**
 * Two-sample ECDF maximum-distance statistic for discrete PAS coordinates.
 * A tie-aware normal approximation is used for the reported p-value.
 */
public final class TieAdjustedKs {
    private TieAdjustedKs() {
    }

    public static final class Result {
        public final double d;
        public final double signedD;
        public final int maxDiffAt0;
        public final double z;
        public final double pvalue;

        Result(
                double distance,
                double signedDistance,
                int maxDifferencePosition0,
                double zScore,
                double pValue
        ) {
            this.d = distance;
            this.signedD = signedDistance;
            this.maxDiffAt0 = maxDifferencePosition0;
            this.z = zScore;
            this.pvalue = pValue;
        }
    }

    /** Counts are keyed by genomic 0-based PAS coordinate. */
    public static Result test(
            Map<Integer, Integer> controlCounts,
            Map<Integer, Integer> caseCounts
    ) {
        TreeSet<Integer> positions = new TreeSet<Integer>();
        positions.addAll(controlCounts.keySet());
        positions.addAll(caseCounts.keySet());

        long controlTotal = sumPositiveCounts(controlCounts);
        long caseTotal = sumPositiveCounts(caseCounts);
        long combinedTotal = controlTotal + caseTotal;

        if (controlTotal <= 0 || caseTotal <= 0 || positions.isEmpty()) {
            return new Result(0.0, 0.0, -1, 0.0, 1.0);
        }

        long cumulativeControl = 0;
        long cumulativeCase = 0;
        double maximumAbsoluteDifference = -1.0;
        double signedDifferenceAtMaximum = 0.0;
        int maximumDifferencePosition = -1;
        double tieVariance = 0.0;

        for (int position : positions) {
            int controlAtPosition = get(controlCounts, position);
            int caseAtPosition = get(caseCounts, position);
            cumulativeControl += controlAtPosition;
            cumulativeCase += caseAtPosition;

            double ecdfDifference =
                    ((double) cumulativeControl / controlTotal)
                            - ((double) cumulativeCase / caseTotal);
            double absoluteDifference = Math.abs(ecdfDifference);

            if (absoluteDifference > maximumAbsoluteDifference) {
                maximumAbsoluteDifference = absoluteDifference;
                signedDifferenceAtMaximum = ecdfDifference;
                maximumDifferencePosition = position;
            }

            double pooledFrequency =
                    (double) (controlAtPosition + caseAtPosition) / combinedTotal;
            tieVariance += pooledFrequency * (1.0 - pooledFrequency);
        }

        double effectiveSampleSize =
                (double) controlTotal * caseTotal / (controlTotal + caseTotal);
        double zScore = 0.0;
        double pValue = 1.0;

        if (tieVariance > 0.0 && effectiveSampleSize > 0.0) {
            zScore = signedDifferenceAtMaximum
                    / Math.sqrt(tieVariance / effectiveSampleSize);
            pValue = MathUtil.erfc(Math.abs(zScore) / Math.sqrt(2.0));
            if (Double.isNaN(pValue) || pValue > 1.0) {
                pValue = 1.0;
            }
            if (pValue < 0.0) {
                pValue = 0.0;
            }
        }

        return new Result(
                Math.max(0.0, maximumAbsoluteDifference),
                signedDifferenceAtMaximum,
                maximumDifferencePosition,
                zScore,
                pValue
        );
    }

    private static int get(Map<Integer, Integer> counts, int position) {
        Integer value = counts.get(position);
        return value == null ? 0 : value;
    }

    private static long sumPositiveCounts(Map<Integer, Integer> counts) {
        long total = 0;
        for (Integer value : counts.values()) {
            if (value != null && value > 0) {
                total += value;
            }
        }
        return total;
    }
}
