package org.apaura.stats;

/** Fisher exact tests for 2x2 contingency tables. */
public final class FisherExact {
    private FisherExact() {
    }

    private static double logHypergeometric(int cellA, int rowOne, int columnOne, int total) {
        return MathUtil.logChoose(columnOne, cellA)
                + MathUtil.logChoose(total - columnOne, rowOne - cellA)
                - MathUtil.logChoose(total, rowOne);
    }

    public static double twoSided(int a, int b, int c, int d) {
        int rowOne = a + b;
        int columnOne = a + c;
        int total = a + b + c + d;
        int lower = Math.max(0, rowOne - (total - columnOne));
        int upper = Math.min(rowOne, columnOne);

        double observedLogProbability = logHypergeometric(a, rowOne, columnOne, total);
        double pValue = 0.0;

        for (int candidateA = lower; candidateA <= upper; candidateA++) {
            double candidateLogProbability =
                    logHypergeometric(candidateA, rowOne, columnOne, total);
            if (candidateLogProbability <= observedLogProbability + 1.0E-10) {
                pValue += Math.exp(candidateLogProbability);
            }
        }
        return Math.min(1.0, pValue);
    }

    public static double rightTail(int a, int b, int c, int d) {
        int rowOne = a + b;
        int columnOne = a + c;
        int total = a + b + c + d;
        int upper = Math.min(rowOne, columnOne);

        double pValue = 0.0;
        for (int candidateA = a; candidateA <= upper; candidateA++) {
            pValue += Math.exp(logHypergeometric(candidateA, rowOne, columnOne, total));
        }
        return Math.min(1.0, pValue);
    }
}
