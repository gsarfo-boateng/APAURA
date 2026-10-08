package org.apaura.rbp;

/** Position-weight matrix parsed from MEME format. */
public final class MemeMotif {
    private static final double MIN_PROBABILITY = 1.0E-9;
    private static final double BACKGROUND_PROBABILITY = 0.25;
    private static final double LOG_TWO = Math.log(2.0);

    public final String id;
    public final String name;
    public final double[][] pwm;
    public final int width;
    public final double minScore;
    public final double maxScore;

    public MemeMotif(String id, String name, double[][] probabilities) {
        this.id = id;
        this.name = name == null || name.isEmpty() ? id : name;
        this.pwm = probabilities;
        this.width = probabilities.length;

        double minimum = 0.0;
        double maximum = 0.0;
        for (double[] row : probabilities) {
            double rowMinimum = Double.POSITIVE_INFINITY;
            double rowMaximum = Double.NEGATIVE_INFINITY;

            for (int baseIndex = 0; baseIndex < 4; baseIndex++) {
                double probability = Math.max(MIN_PROBABILITY, row[baseIndex]);
                double logOdds = log2(probability / BACKGROUND_PROBABILITY);
                rowMinimum = Math.min(rowMinimum, logOdds);
                rowMaximum = Math.max(rowMaximum, logOdds);
            }

            minimum += rowMinimum;
            maximum += rowMaximum;
        }

        this.minScore = minimum;
        this.maxScore = maximum;
    }

    public double normalizedScore(String sequence, int start) {
        double score = 0.0;
        for (int position = 0; position < width; position++) {
            int baseIndex = baseIndex(sequence.charAt(start + position));
            if (baseIndex < 0) {
                return Double.NEGATIVE_INFINITY;
            }

            double probability = Math.max(
                    MIN_PROBABILITY,
                    pwm[position][baseIndex]
            );
            score += log2(probability / BACKGROUND_PROBABILITY);
        }

        if (maxScore <= minScore) {
            return 0.0;
        }
        return (score - minScore) / (maxScore - minScore);
    }

    private static double log2(double value) {
        return Math.log(value) / LOG_TWO;
    }

    private static int baseIndex(char base) {
        switch (Character.toUpperCase(base)) {
            case 'A':
                return 0;
            case 'C':
                return 1;
            case 'G':
                return 2;
            case 'T':
            case 'U':
                return 3;
            default:
                return -1;
        }
    }
}
