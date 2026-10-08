package org.apaura.stats;

/** Weighted pool-adjacent-violators regression constrained to y0 >= y1 >= ... >= 0. */
public final class IsotonicRegression {
    private IsotonicRegression() {
    }

    public static double[] decreasing(double[] values, double[] weights) {
        int length = values.length;
        if (weights.length != length) {
            throw new IllegalArgumentException("weight length mismatch");
        }

        double[] blockMeans = new double[length];
        double[] blockWeights = new double[length];
        int[] blockStarts = new int[length];
        int[] blockEnds = new int[length];
        int blockCount = 0;

        for (int i = 0; i < length; i++) {
            blockMeans[blockCount] = Math.max(0.0, values[i]);
            blockWeights[blockCount] = Math.max(1.0E-12, weights[i]);
            blockStarts[blockCount] = i;
            blockEnds[blockCount] = i;
            blockCount++;

            while (blockCount >= 2
                    && blockMeans[blockCount - 2] < blockMeans[blockCount - 1]) {
                double mergedWeight =
                        blockWeights[blockCount - 2] + blockWeights[blockCount - 1];
                double mergedMean =
                        (blockMeans[blockCount - 2] * blockWeights[blockCount - 2]
                                + blockMeans[blockCount - 1] * blockWeights[blockCount - 1])
                                / mergedWeight;

                blockMeans[blockCount - 2] = mergedMean;
                blockWeights[blockCount - 2] = mergedWeight;
                blockEnds[blockCount - 2] = blockEnds[blockCount - 1];
                blockCount--;
            }
        }

        double[] fitted = new double[length];
        for (int block = 0; block < blockCount; block++) {
            for (int i = blockStarts[block]; i <= blockEnds[block]; i++) {
                fitted[i] = Math.max(0.0, blockMeans[block]);
            }
        }
        return fitted;
    }
}
