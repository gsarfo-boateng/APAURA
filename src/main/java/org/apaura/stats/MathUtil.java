package org.apaura.stats;

/** Mathematical helper functions used by APAURA statistical tests. */
public final class MathUtil {
    private static final double[] LANCZOS_COEFFICIENTS = {
            676.5203681218851,
            -1259.1392167224028,
            771.3234287776531,
            -176.6150291621406,
            12.507343278686905,
            -0.13857109526572012,
            9.984369578019572E-6,
            1.5056327351493116E-7
    };

    private MathUtil() {
    }

    public static double logGamma(double z) {
        if (z < 0.5) {
            return Math.log(Math.PI)
                    - Math.log(Math.sin(Math.PI * z))
                    - logGamma(1.0 - z);
        }

        double shifted = z - 1.0;
        double series = 0.9999999999998099;
        for (int i = 0; i < LANCZOS_COEFFICIENTS.length; i++) {
            series += LANCZOS_COEFFICIENTS[i] / (shifted + i + 1.0);
        }

        double t = shifted + LANCZOS_COEFFICIENTS.length - 0.5;
        return 0.5 * Math.log(2.0 * Math.PI)
                + (shifted + 0.5) * Math.log(t)
                - t
                + Math.log(series);
    }

    public static double logChoose(int n, int k) {
        if (k < 0 || k > n) {
            return Double.NEGATIVE_INFINITY;
        }
        return logGamma(n + 1.0)
                - logGamma(k + 1.0)
                - logGamma(n - k + 1.0);
    }

    public static double erfc(double x) {
        double absoluteX = Math.abs(x);
        double t = 1.0 / (1.0 + 0.5 * absoluteX);
        double approximation = t * Math.exp(
                -absoluteX * absoluteX
                        - 1.26551223
                        + t * (1.00002368
                        + t * (0.37409196
                        + t * (0.09678418
                        + t * (-0.18628806
                        + t * (0.27886807
                        + t * (-1.13520398
                        + t * (1.48851587
                        + t * (-0.82215223
                        + t * 0.17087277))))))))
        );
        return x >= 0.0 ? approximation : 2.0 - approximation;
    }

    public static double chiSquare1Survival(double statistic) {
        if (statistic <= 0.0) {
            return 1.0;
        }
        return erfc(Math.sqrt(statistic / 2.0));
    }
}
