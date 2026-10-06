package org.apaura.stats;

import java.util.Arrays;
import java.util.Comparator;

/** Multiple-testing correction utilities. */
public final class MultipleTesting {
    private MultipleTesting() {
    }

    /** Benjamini-Hochberg false-discovery-rate correction. */
    public static double[] bh(final double[] pValues) {
        Integer[] order = new Integer[pValues.length];
        for (int i = 0; i < pValues.length; i++) {
            order[i] = i;
        }

        Arrays.sort(order, new Comparator<Integer>() {
            @Override
            public int compare(Integer left, Integer right) {
                return Double.compare(pValues[left], pValues[right]);
            }
        });

        double[] qValues = new double[pValues.length];
        double runningMinimum = 1.0;

        for (int rankIndex = pValues.length - 1; rankIndex >= 0; rankIndex--) {
            int originalIndex = order[rankIndex];
            double adjusted = pValues[originalIndex] * pValues.length / (rankIndex + 1.0);
            runningMinimum = Math.min(runningMinimum, adjusted);
            qValues[originalIndex] = Math.min(1.0, runningMinimum);
        }

        return qValues;
    }
}
