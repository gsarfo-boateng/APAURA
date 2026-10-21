package org.apaura.stats;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Fixed-precision Dirichlet-multinomial likelihood engine.
 *
 * <p>This class is deliberately independent of precision estimation/moderation.
 * It provides the numerical core that the APAURA moderated-DM layer will use.
 * It is not wired into production short-read inference until moderated
 * precision estimation has been implemented and validated against the frozen
 * DRIMSeq reference.</p>
 */
public final class DirichletMultinomial {

    private static final double MIN_PROP = 1.0E-12;
    private static final double GRAD_TOL = 1.0E-8;
    private static final double STEP_TOL = 1.0E-12;
    private static final int MAX_ITER = 500;

    private DirichletMultinomial() {
        // Utility class.
    }

    public static final class FitResult {
        public final double[] proportions;
        public final double logLikelihood;
        public final boolean converged;
        public final int iterations;

        FitResult(
                double[] proportions,
                double logLikelihood,
                boolean converged,
                int iterations
        ) {
            this.proportions = proportions;
            this.logLikelihood = logLikelihood;
            this.converged = converged;
            this.iterations = iterations;
        }
    }

    public static final class TestResult {
        public final double precision;
        public final double nullLogLikelihood;
        public final double fullLogLikelihood;
        public final double likelihoodRatio;
        public final int degreesOfFreedom;
        public final double pvalue;
        public final boolean converged;

        TestResult(
                double precision,
                double nullLogLikelihood,
                double fullLogLikelihood,
                double likelihoodRatio,
                int degreesOfFreedom,
                double pvalue,
                boolean converged
        ) {
            this.precision = precision;
            this.nullLogLikelihood = nullLogLikelihood;
            this.fullLogLikelihood = fullLogLikelihood;
            this.likelihoodRatio = likelihoodRatio;
            this.degreesOfFreedom = degreesOfFreedom;
            this.pvalue = pvalue;
            this.converged = converged;
        }
    }

    /**
     * Fit one composition vector for a sample-by-feature count matrix at a
     * fixed precision.
     */
    public static FitResult fitOneGroup(int[][] counts, double precision) {
        validateCounts(counts);
        if (!(precision > 0.0) || !Double.isFinite(precision)) {
            return failedFit(counts[0].length);
        }

        final int originalFeatures = counts[0].length;

        List<Integer> nonzeroSamples = new ArrayList<Integer>();
        for (int sample = 0; sample < counts.length; sample++) {
            if (rowSum(counts[sample]) > 0) {
                nonzeroSamples.add(sample);
            }
        }
        if (nonzeroSamples.isEmpty()) {
            return failedFit(originalFeatures);
        }

        long[] featureTotals = new long[originalFeatures];
        for (Integer sample : nonzeroSamples) {
            for (int feature = 0; feature < originalFeatures; feature++) {
                featureTotals[feature] += counts[sample][feature];
            }
        }

        int activeCount = 0;
        for (long total : featureTotals) {
            if (total > 0L) {
                activeCount++;
            }
        }

        // At least two active competing features are required. The reduced
        // matrix below chooses its own active reference feature, so a zero
        // count in the original final column does not make the event invalid.
        if (activeCount < 2) {
            return failedFit(originalFeatures);
        }

        int[] active = new int[activeCount];
        int cursor = 0;
        for (int feature = 0; feature < originalFeatures; feature++) {
            if (featureTotals[feature] > 0L) {
                active[cursor++] = feature;
            }
        }

        int[][] reduced = new int[nonzeroSamples.size()][activeCount];
        for (int i = 0; i < nonzeroSamples.size(); i++) {
            int sample = nonzeroSamples.get(i);
            for (int j = 0; j < activeCount; j++) {
                reduced[i][j] = counts[sample][active[j]];
            }
        }

        double[] initial = empiricalProportions(reduced);

        if (reduced.length == 1) {
            double ll = logLikelihood(reduced, initial, precision);
            return expandFit(initial, active, originalFeatures, ll, true, 0);
        }

        final int dim = activeCount - 1;
        double[] eta = new double[dim];
        double reference = Math.max(initial[activeCount - 1], MIN_PROP);
        for (int j = 0; j < dim; j++) {
            eta[j] = Math.log(Math.max(initial[j], MIN_PROP) / reference);
        }

        double[][] inverseHessian = identity(dim);
        Objective current = objective(reduced, precision, eta);
        if (!Double.isFinite(current.value) || !allFinite(current.gradient)) {
            return failedFit(originalFeatures);
        }

        boolean converged = false;
        int iterations = 0;

        for (iterations = 0; iterations < MAX_ITER; iterations++) {
            if (maxAbs(current.gradient) < GRAD_TOL) {
                converged = true;
                break;
            }

            double[] direction = matVec(inverseHessian, current.gradient);
            double directionalDerivative = dot(current.gradient, direction);

            // A positive directional derivative is required for maximization.
            // If BFGS loses positive-definiteness, restart with steepest ascent.
            if (!(directionalDerivative > 0.0) || !Double.isFinite(directionalDerivative)) {
                inverseHessian = identity(dim);
                direction = Arrays.copyOf(current.gradient, dim);
                directionalDerivative = dot(current.gradient, direction);
            }

            double step = 1.0;
            Objective candidate = null;
            double[] etaCandidate = null;
            final double armijo = 1.0E-4;

            while (step >= STEP_TOL) {
                etaCandidate = addScaled(eta, direction, step);
                candidate = objective(reduced, precision, etaCandidate);

                if (Double.isFinite(candidate.value)
                        && allFinite(candidate.gradient)
                        && candidate.value >= current.value
                        + armijo * step * directionalDerivative) {
                    break;
                }
                step *= 0.5;
            }

            if (step < STEP_TOL || candidate == null || etaCandidate == null) {
                break;
            }

            double[] s = subtract(etaCandidate, eta);
            double[] y = subtract(candidate.gradient, current.gradient);

            // Standard BFGS is stated for minimization. For maximization,
            // update using y = oldGradient - newGradient so that the local
            // curvature matrix represents -H.
            double[] yForMax = scale(y, -1.0);
            double sy = dot(s, yForMax);

            if (sy > 1.0E-12 && Double.isFinite(sy)) {
                inverseHessian = bfgsInverseUpdate(inverseHessian, s, yForMax, sy);
            } else {
                inverseHessian = identity(dim);
            }

            eta = etaCandidate;
            current = candidate;

            if (maxAbs(s) < STEP_TOL) {
                converged = true;
                break;
            }
        }

        if (!converged && maxAbs(current.gradient) < 1.0E-6) {
            converged = true;
        }

        double[] fitted = softmaxReference(eta);
        return expandFit(
                fitted,
                active,
                originalFeatures,
                current.value,
                converged,
                iterations
        );
    }

    /**
     * Full-vs-null two-condition DM likelihood-ratio test at a fixed precision.
     *
     * <p>groups must contain only 0, 1, or -1. Samples marked -1 are omitted.</p>
     */
    public static TestResult testTwoGroup(
            int[][] counts,
            int[] groups,
            double precision
    ) {
        validateCounts(counts);
        if (groups == null || groups.length != counts.length) {
            throw new IllegalArgumentException(
                    "groups length must equal number of samples"
            );
        }

        List<int[]> all = new ArrayList<int[]>();
        List<int[]> group0 = new ArrayList<int[]>();
        List<int[]> group1 = new ArrayList<int[]>();

        for (int i = 0; i < counts.length; i++) {
            if (groups[i] == -1) {
                continue;
            }
            if (groups[i] != 0 && groups[i] != 1) {
                throw new IllegalArgumentException(
                        "groups must contain only -1, 0, or 1"
                );
            }
            if (rowSum(counts[i]) <= 0) {
                continue;
            }
            int[] row = Arrays.copyOf(counts[i], counts[i].length);
            all.add(row);
            if (groups[i] == 0) {
                group0.add(row);
            } else {
                group1.add(row);
            }
        }

        int features = counts[0].length;
        int df = Math.max(0, features - 1);

        if (group0.isEmpty() || group1.isEmpty() || df == 0) {
            return failedTest(precision, df);
        }

        FitResult nullFit = fitOneGroup(toMatrix(all), precision);
        FitResult fit0 = fitOneGroup(toMatrix(group0), precision);
        FitResult fit1 = fitOneGroup(toMatrix(group1), precision);

        /*
         * Use the best finite likelihood returned by the optimizer even when
         * its strict gradient/step convergence flag is false. Numerical
         * validity is determined by the likelihood/LR/p-value checks below.
         */
        if (!Double.isFinite(nullFit.logLikelihood)
                || !Double.isFinite(fit0.logLikelihood)
                || !Double.isFinite(fit1.logLikelihood)) {
            return failedTest(precision, df);
        }

        double full = fit0.logLikelihood + fit1.logLikelihood;
        double lr = 2.0 * (full - nullFit.logLikelihood);

        // Numerical optimization can create tiny negative LRs under the null.
        if (lr < 0.0 && lr > -1.0E-8) {
            lr = 0.0;
        }
        if (lr < 0.0 || !Double.isFinite(lr)) {
            return failedTest(precision, df);
        }

        double p = chiSquareSurvival(lr, df);

        return new TestResult(
                precision,
                nullFit.logLikelihood,
                full,
                lr,
                df,
                p,
                Double.isFinite(p)
        );
    }

    /**
     * DM log-likelihood without the multinomial normalizing constant. The
     * omitted term is independent of proportions and precision and cancels in
     * likelihood-ratio tests.
     */
    public static double logLikelihood(
            int[][] counts,
            double[] proportions,
            double precision
    ) {
        validateCounts(counts);
        if (proportions == null || proportions.length != counts[0].length) {
            throw new IllegalArgumentException(
                    "proportions length must equal number of features"
            );
        }
        if (!(precision > 0.0) || !Double.isFinite(precision)) {
            return Double.NaN;
        }

        double sumP = 0.0;
        for (double p : proportions) {
            if (!(p > 0.0) || !Double.isFinite(p)) {
                return Double.NaN;
            }
            sumP += p;
        }
        if (Math.abs(sumP - 1.0) > 1.0E-8) {
            return Double.NaN;
        }

        double ll = 0.0;
        for (int[] sample : counts) {
            int total = rowSum(sample);
            if (total <= 0) {
                continue;
            }
            ll += logGamma(precision) - logGamma(total + precision);
            for (int j = 0; j < sample.length; j++) {
                double alpha = precision * proportions[j];
                ll += logGamma(sample[j] + alpha) - logGamma(alpha);
            }
        }
        return ll;
    }

    private static Objective objective(
            int[][] counts,
            double precision,
            double[] eta
    ) {
        double[] p = softmaxReference(eta);
        double value = logLikelihood(counts, p, precision);

        double[] independentGradient = new double[p.length];
        for (int feature = 0; feature < p.length; feature++) {
            double alpha = precision * p[feature];
            double g = 0.0;
            for (int[] sample : counts) {
                if (rowSum(sample) <= 0) {
                    continue;
                }
                g += digamma(sample[feature] + alpha) - digamma(alpha);
            }
            independentGradient[feature] = precision * g;
        }

        double weighted = 0.0;
        for (int feature = 0; feature < p.length; feature++) {
            weighted += p[feature] * independentGradient[feature];
        }

        double[] gradient = new double[eta.length];
        for (int feature = 0; feature < eta.length; feature++) {
            gradient[feature] =
                    p[feature] * (independentGradient[feature] - weighted);
        }

        return new Objective(value, gradient);
    }

    private static double[] empiricalProportions(int[][] counts) {
        long[] totals = new long[counts[0].length];
        long grand = 0L;
        for (int[] sample : counts) {
            for (int j = 0; j < sample.length; j++) {
                totals[j] += sample[j];
                grand += sample[j];
            }
        }

        double[] p = new double[totals.length];
        double denominator = grand + MIN_PROP * totals.length;
        for (int j = 0; j < totals.length; j++) {
            p[j] = (totals[j] + MIN_PROP) / denominator;
        }
        renormalize(p);
        return p;
    }

    private static FitResult expandFit(
            double[] activeProportions,
            int[] active,
            int originalFeatures,
            double ll,
            boolean converged,
            int iterations
    ) {
        double[] expanded = new double[originalFeatures];
        for (int j = 0; j < active.length; j++) {
            expanded[active[j]] = activeProportions[j];
        }
        return new FitResult(expanded, ll, converged, iterations);
    }

    private static FitResult failedFit(int features) {
        double[] p = new double[features];
        Arrays.fill(p, Double.NaN);
        return new FitResult(p, Double.NaN, false, 0);
    }

    private static TestResult failedTest(double precision, int df) {
        return new TestResult(
                precision,
                Double.NaN,
                Double.NaN,
                Double.NaN,
                df,
                Double.NaN,
                false
        );
    }

    private static double[] softmaxReference(double[] eta) {
        double max = 0.0; // reference logit is exactly zero
        for (double v : eta) {
            if (v > max) {
                max = v;
            }
        }

        double denominator = Math.exp(-max);
        double[] p = new double[eta.length + 1];

        for (int j = 0; j < eta.length; j++) {
            p[j] = Math.exp(eta[j] - max);
            denominator += p[j];
        }

        for (int j = 0; j < eta.length; j++) {
            p[j] /= denominator;
            p[j] = Math.max(p[j], MIN_PROP);
        }
        p[p.length - 1] = Math.max(Math.exp(-max) / denominator, MIN_PROP);
        renormalize(p);
        return p;
    }

    private static void renormalize(double[] p) {
        double sum = 0.0;
        for (double v : p) {
            sum += v;
        }
        for (int j = 0; j < p.length; j++) {
            p[j] /= sum;
        }
    }

    private static void validateCounts(int[][] counts) {
        if (counts == null || counts.length == 0
                || counts[0] == null || counts[0].length < 2) {
            throw new IllegalArgumentException(
                    "counts must contain at least one sample and two features"
            );
        }

        int features = counts[0].length;
        for (int i = 0; i < counts.length; i++) {
            if (counts[i] == null || counts[i].length != features) {
                throw new IllegalArgumentException(
                        "counts must be a rectangular sample-by-feature matrix"
                );
            }
            for (int count : counts[i]) {
                if (count < 0) {
                    throw new IllegalArgumentException(
                            "counts must be non-negative"
                    );
                }
            }
        }
    }

    private static int rowSum(int[] row) {
        int total = 0;
        for (int value : row) {
            total += value;
        }
        return total;
    }

    private static int[][] toMatrix(List<int[]> rows) {
        int[][] out = new int[rows.size()][];
        for (int i = 0; i < rows.size(); i++) {
            out[i] = rows.get(i);
        }
        return out;
    }

    private static final class Objective {
        final double value;
        final double[] gradient;

        Objective(double value, double[] gradient) {
            this.value = value;
            this.gradient = gradient;
        }
    }

    private static double[][] identity(int n) {
        double[][] out = new double[n][n];
        for (int i = 0; i < n; i++) {
            out[i][i] = 1.0;
        }
        return out;
    }

    private static double[] matVec(double[][] matrix, double[] vector) {
        double[] out = new double[vector.length];
        for (int i = 0; i < matrix.length; i++) {
            double value = 0.0;
            for (int j = 0; j < vector.length; j++) {
                value += matrix[i][j] * vector[j];
            }
            out[i] = value;
        }
        return out;
    }

    private static double[][] bfgsInverseUpdate(
            double[][] h,
            double[] s,
            double[] y,
            double sy
    ) {
        int n = s.length;
        double rho = 1.0 / sy;
        double[][] left = identity(n);
        double[][] right = identity(n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                left[i][j] -= rho * s[i] * y[j];
                right[i][j] -= rho * y[i] * s[j];
            }
        }

        double[][] temp = multiply(left, h);
        double[][] out = multiply(temp, right);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                out[i][j] += rho * s[i] * s[j];
            }
        }
        return out;
    }

    private static double[][] multiply(double[][] a, double[][] b) {
        int n = a.length;
        double[][] out = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < n; k++) {
                double aik = a[i][k];
                for (int j = 0; j < n; j++) {
                    out[i][j] += aik * b[k][j];
                }
            }
        }
        return out;
    }

    private static double[] addScaled(double[] a, double[] b, double scale) {
        double[] out = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            out[i] = a[i] + scale * b[i];
        }
        return out;
    }

    private static double[] subtract(double[] a, double[] b) {
        double[] out = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            out[i] = a[i] - b[i];
        }
        return out;
    }

    private static double[] scale(double[] a, double scale) {
        double[] out = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            out[i] = scale * a[i];
        }
        return out;
    }

    private static double dot(double[] a, double[] b) {
        double out = 0.0;
        for (int i = 0; i < a.length; i++) {
            out += a[i] * b[i];
        }
        return out;
    }

    private static double maxAbs(double[] values) {
        double out = 0.0;
        for (double value : values) {
            out = Math.max(out, Math.abs(value));
        }
        return out;
    }

    private static boolean allFinite(double[] values) {
        for (double value : values) {
            if (!Double.isFinite(value)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Lanczos log-gamma approximation with reflection for x < 0.5.
     */
    static double logGamma(double x) {
        final double[] coefficients = {
                0.99999999999980993,
                676.5203681218851,
                -1259.1392167224028,
                771.32342877765313,
                -176.61502916214059,
                12.507343278686905,
                -0.13857109526572012,
                9.9843695780195716E-6,
                1.5056327351493116E-7
        };

        if (x < 0.5) {
            return Math.log(Math.PI)
                    - Math.log(Math.sin(Math.PI * x))
                    - logGamma(1.0 - x);
        }

        double z = x - 1.0;
        double a = coefficients[0];
        for (int i = 1; i < coefficients.length; i++) {
            a += coefficients[i] / (z + i);
        }

        double t = z + 7.5;
        return 0.5 * Math.log(2.0 * Math.PI)
                + (z + 0.5) * Math.log(t)
                - t
                + Math.log(a);
    }

    /**
     * Digamma approximation accurate for positive arguments used in the DM
     * score function.
     */
    static double digamma(double x) {
        if (!(x > 0.0) || !Double.isFinite(x)) {
            return Double.NaN;
        }

        double result = 0.0;
        while (x < 8.0) {
            result -= 1.0 / x;
            x += 1.0;
        }

        double inv = 1.0 / x;
        double inv2 = inv * inv;

        result += Math.log(x)
                - 0.5 * inv
                - inv2 * (
                        1.0 / 12.0
                        - inv2 * (
                                1.0 / 120.0
                                - inv2 * (
                                        1.0 / 252.0
                                        - inv2 * (
                                                1.0 / 240.0
                                                - inv2 * 5.0 / 660.0
                                        )
                                )
                        )
                );
        return result;
    }

    static double chiSquareSurvival(double statistic, int df) {
        if (statistic < 0.0 || df <= 0 || !Double.isFinite(statistic)) {
            return Double.NaN;
        }
        if (statistic == 0.0) {
            return 1.0;
        }
        return regularizedGammaQ(0.5 * df, 0.5 * statistic);
    }

    /**
     * Regularized upper incomplete gamma Q(a,x).
     */
    static double regularizedGammaQ(double a, double x) {
        if (!(a > 0.0) || x < 0.0 || Double.isNaN(x)) {
            return Double.NaN;
        }
        if (x == 0.0) {
            return 1.0;
        }

        if (x < a + 1.0) {
            return clamp01(1.0 - regularizedGammaPSeries(a, x));
        }

        final double tiny = 1.0E-300;
        final double epsilon = 1.0E-14;
        final int maxIterations = 10000;

        double b = x + 1.0 - a;
        double c = 1.0 / tiny;
        double d = 1.0 / Math.max(Math.abs(b), tiny);
        if (b < 0.0) {
            d = -d;
        }
        double h = d;

        for (int i = 1; i <= maxIterations; i++) {
            double an = -i * (i - a);
            b += 2.0;

            d = an * d + b;
            if (Math.abs(d) < tiny) {
                d = tiny;
            }

            c = b + an / c;
            if (Math.abs(c) < tiny) {
                c = tiny;
            }

            d = 1.0 / d;
            double delta = d * c;
            h *= delta;

            if (Math.abs(delta - 1.0) < epsilon) {
                break;
            }
        }

        double front = Math.exp(
                -x + a * Math.log(x) - logGamma(a)
        );
        return clamp01(front * h);
    }

    private static double regularizedGammaPSeries(double a, double x) {
        final double epsilon = 1.0E-14;
        final int maxIterations = 10000;

        double sum = 1.0 / a;
        double term = sum;
        double ap = a;

        for (int n = 1; n <= maxIterations; n++) {
            ap += 1.0;
            term *= x / ap;
            sum += term;
            if (Math.abs(term) < Math.abs(sum) * epsilon) {
                break;
            }
        }

        return clamp01(
                sum * Math.exp(-x + a * Math.log(x) - logGamma(a))
        );
    }

    private static double clamp01(double value) {
        if (value < 0.0) {
            return 0.0;
        }
        if (value > 1.0) {
            return 1.0;
        }
        return value;
    }
}
