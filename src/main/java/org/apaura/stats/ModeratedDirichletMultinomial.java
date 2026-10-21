package org.apaura.stats;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Replicate-aware, event-level Dirichlet-multinomial inference for APAURA.
 *
 * <p>The implementation follows the statistical structure calibrated against
 * DRIMSeq: Cox-Reid adjusted profile likelihood, a common precision, a
 * 21-point log2 genewise precision grid, expression-trended likelihood
 * moderation, and a full-vs-null DM likelihood-ratio test. It is implemented
 * independently in Java and has no R runtime dependency.</p>
 */
public final class ModeratedDirichletMultinomial {

    private static final double MIN_PRECISION = 1.0E-4;
    private static final double MAX_COMMON_PRECISION = 1000.0;
    private static final double COMMON_PRECISION_TOL = 10.0;
    private static final int GRID_LENGTH = 21;
    private static final int GRID_MIN = -10;
    private static final double TREND_SPAN = 0.10;
    private static final double LOESS_SPAN = 0.75;

    private ModeratedDirichletMultinomial() {
    }

    public static final class Event {
        public final String id;
        public final int[][] counts;
        public final int[] groups;
        public final double meanExpression;

        public Event(String id, int[][] counts, int[] groups) {
            if (id == null || id.isEmpty()) {
                throw new IllegalArgumentException("event id must not be empty");
            }
            if (counts == null || counts.length == 0
                    || groups == null || groups.length != counts.length) {
                throw new IllegalArgumentException(
                        "counts and groups must contain the same non-zero number of samples"
                );
            }
            if (counts[0] == null || counts[0].length < 2) {
                throw new IllegalArgumentException(
                        "DM events require at least two competing features"
                );
            }

            int features = counts[0].length;
            this.counts = new int[counts.length][features];
            this.groups = Arrays.copyOf(groups, groups.length);
            long total = 0L;

            for (int i = 0; i < counts.length; i++) {
                if (counts[i] == null || counts[i].length != features) {
                    throw new IllegalArgumentException("counts must be rectangular");
                }
                if (groups[i] != 0 && groups[i] != 1) {
                    throw new IllegalArgumentException("groups must contain only 0 or 1");
                }
                for (int j = 0; j < features; j++) {
                    if (counts[i][j] < 0) {
                        throw new IllegalArgumentException("counts must be non-negative");
                    }
                    this.counts[i][j] = counts[i][j];
                    total += counts[i][j];
                }
            }

            this.id = id;
            this.meanExpression = total / (double) counts.length;
        }
    }

    public static final class Result {
        public final String id;
        public final double precision;
        public final double likelihoodRatio;
        public final int degreesOfFreedom;
        public final double pvalue;
        public final boolean converged;

        Result(
                String id,
                double precision,
                double likelihoodRatio,
                int degreesOfFreedom,
                double pvalue,
                boolean converged
        ) {
            this.id = id;
            this.precision = precision;
            this.likelihoodRatio = likelihoodRatio;
            this.degreesOfFreedom = degreesOfFreedom;
            this.pvalue = pvalue;
            this.converged = converged;
        }
    }

    public static final class BatchResult {
        public final double commonPrecision;
        public final List<Result> results;
        private final Map<String, Result> byId;

        BatchResult(double commonPrecision, List<Result> results) {
            this.commonPrecision = commonPrecision;
            this.results = Collections.unmodifiableList(results);
            this.byId = new LinkedHashMap<String, Result>();
            for (Result result : results) {
                this.byId.put(result.id, result);
            }
        }

        public Result get(String id) {
            return byId.get(id);
        }
    }

    public static BatchResult test(List<Event> events) {
        if (events == null || events.isEmpty()) {
            return new BatchResult(Double.NaN, new ArrayList<Result>());
        }

        List<Event> usable = new ArrayList<Event>();
        for (Event event : events) {
            if (hasBothGroups(event) && event.counts[0].length >= 2) {
                usable.add(event);
            }
        }
        if (usable.isEmpty()) {
            return new BatchResult(Double.NaN, new ArrayList<Result>());
        }

        double commonPrecision = estimateCommonPrecision(usable);
        if (!(commonPrecision > 0.0) || !Double.isFinite(commonPrecision)) {
            commonPrecision = 100.0;
        }

        double[] precisions = estimateModeratedPrecisions(usable, commonPrecision);

        List<Result> results = new ArrayList<Result>();
        for (int i = 0; i < usable.size(); i++) {
            Event event = usable.get(i);
            double precision = precisions[i];

            /*
             * Match the calibrated DRIMSeq missingness semantics. Genewise
             * precision is defined only for events with a complete 21-point
             * adjusted profile-likelihood grid. DRIMSeq leaves all other
             * events NA; it does not rescue them with the common precision.
             */
            if (!(precision > 0.0) || !Double.isFinite(precision)) {
                results.add(new Result(
                        event.id,
                        Double.NaN,
                        Double.NaN,
                        Math.max(1, event.counts[0].length - 1),
                        Double.NaN,
                        false
                ));
                continue;
            }

            DirichletMultinomial.TestResult test =
                    DirichletMultinomial.testTwoGroup(
                            event.counts,
                            event.groups,
                            precision
                    );

            /*
             * Likewise, do not substitute the common precision if the
             * genewise fixed-precision LRT is numerically invalid.
             */
            boolean valid = Double.isFinite(test.pvalue)
                    && Double.isFinite(test.likelihoodRatio);

            results.add(new Result(
                    event.id,
                    precision,
                    test.likelihoodRatio,
                    test.degreesOfFreedom,
                    test.pvalue,
                    valid
            ));
        }

        return new BatchResult(commonPrecision, results);
    }

    private static boolean hasBothGroups(Event event) {
        boolean zero = false;
        boolean one = false;
        for (int group : event.groups) {
            zero |= group == 0;
            one |= group == 1;
        }
        return zero && one;
    }

    private static double estimateCommonPrecision(List<Event> events) {
        List<Event> candidates = new ArrayList<Event>();
        for (Event event : events) {
            double value = adjustedProfileLikelihood(event, 100.0);
            if (Double.isFinite(value)) {
                candidates.add(event);
            }
        }
        if (candidates.isEmpty()) {
            return Double.NaN;
        }

        if (candidates.size() >= 50) {
            Collections.sort(candidates, new Comparator<Event>() {
                @Override
                public int compare(Event left, Event right) {
                    int c = Long.compareUnsigned(
                            stableHash(left.id),
                            stableHash(right.id)
                    );
                    return c != 0 ? c : left.id.compareTo(right.id);
                }
            });
            int keep = Math.max(
                    1,
                    (int) Math.ceil(candidates.size() * 0.10)
            );
            candidates = new ArrayList<Event>(
                    candidates.subList(0, keep)
            );
        }

        final List<Event> subset = candidates;
        double a = MIN_PRECISION;
        double b = MAX_COMMON_PRECISION;
        final double ratio = 0.5 * (Math.sqrt(5.0) - 1.0);
        double c = b - ratio * (b - a);
        double d = a + ratio * (b - a);
        double fc = commonObjective(subset, c);
        double fd = commonObjective(subset, d);

        for (int iter = 0; iter < 80 && (b - a) > COMMON_PRECISION_TOL; iter++) {
            if (fc > fd) {
                b = d;
                d = c;
                fd = fc;
                c = b - ratio * (b - a);
                fc = commonObjective(subset, c);
            } else {
                a = c;
                c = d;
                fc = fd;
                d = a + ratio * (b - a);
                fd = commonObjective(subset, d);
            }
        }

        double out = 0.5 * (a + b);
        return Double.isFinite(out) && out > 0.0 ? out : Double.NaN;
    }

    private static double commonObjective(List<Event> events, double precision) {
        double sum = 0.0;
        int finite = 0;
        for (Event event : events) {
            double value = adjustedProfileLikelihood(event, precision);
            if (Double.isFinite(value)) {
                sum += value;
                finite++;
            }
        }
        if (finite == 0) {
            return Double.NEGATIVE_INFINITY;
        }
        return (sum / finite) * events.size();
    }

    private static double[] estimateModeratedPrecisions(
            List<Event> events,
            double commonPrecision
    ) {
        int n = events.size();
        double[][] logLik = new double[n][GRID_LENGTH];
        boolean[] complete = new boolean[n];

        for (int i = 0; i < n; i++) {
            complete[i] = true;
            for (int j = 0; j < GRID_LENGTH; j++) {
                double precision = commonPrecision
                        * Math.pow(2.0, GRID_MIN + j);
                logLik[i][j] = adjustedProfileLikelihood(
                        events.get(i),
                        precision
                );
                if (!Double.isFinite(logLik[i][j])) {
                    complete[i] = false;
                }
            }
        }

        List<Integer> valid = new ArrayList<Integer>();
        for (int i = 0; i < n; i++) {
            if (complete[i]) valid.add(i);
        }

        double[] out = new double[n];
        Arrays.fill(out, Double.NaN);
        if (valid.isEmpty()) {
            return out;
        }

        int m = valid.size();
        double[][] table = new double[m][GRID_LENGTH];
        double[] expression = new double[m];

        for (int r = 0; r < m; r++) {
            int source = valid.get(r);
            table[r] = Arrays.copyOf(logLik[source], GRID_LENGTH);
            expression[r] = events.get(source).meanExpression;
        }

        moderateProfileLikelihood(table, expression);

        for (int r = 0; r < m; r++) {
            double x = maximizeGrid(table[r]);
            double precision = commonPrecision * Math.pow(2.0, x);
            if (Double.isFinite(precision) && precision > 0.0) {
                out[valid.get(r)] = precision;
            }
        }

        return out;
    }

    private static void moderateProfileLikelihood(
            double[][] logLik,
            double[] expression
    ) {
        int n = logLik.length;
        if (n == 0) return;

        int[] maxIndex = new int[n];
        boolean[] interior = new boolean[n];
        boolean[] upperBoundary = new boolean[n];
        int interiorCount = 0;
        int upperCount = 0;

        for (int i = 0; i < n; i++) {
            maxIndex[i] = argMax(logLik[i]);
            interior[i] = maxIndex[i] > 0
                    && maxIndex[i] < GRID_LENGTH - 4;
            upperBoundary[i] = maxIndex[i] == GRID_LENGTH - 1;
            if (interior[i]) interiorCount++;
            if (upperBoundary[i]) upperCount++;
        }

        if (interiorCount == 0) return;

        double[][] trend = movingAverageTrend(
                logLik,
                expression,
                interior
        );

        double[] priorN = new double[n];
        if (upperCount > 10) {
            List<Double> boundaryX = new ArrayList<Double>();
            List<Double> boundaryY = new ArrayList<Double>();
            List<Double> scalarPrior = new ArrayList<Double>();

            for (int i = 0; i < n; i++) {
                if (!upperBoundary[i]) continue;

                double eventSpan = span(logLik[i]);
                double trendSpan = span(trend[i]);
                if (!(eventSpan > 0.0)
                        || !(trendSpan > 0.0)
                        || !(expression[i] > 0.0)) {
                    continue;
                }

                double value = eventSpan / trendSpan;
                if (!(value > 0.0) || !Double.isFinite(value)) {
                    continue;
                }

                scalarPrior.add(value);
                boundaryX.add(Math.log10(expression[i]));
                boundaryY.add(Math.log10(value));
            }

            if (scalarPrior.size() > 100) {
                double fallback = median(scalarPrior);
                for (int i = 0; i < n; i++) {
                    double predicted = localQuadraticLoess(
                            boundaryX,
                            boundaryY,
                            Math.log10(Math.max(expression[i], 1.0E-12)),
                            LOESS_SPAN
                    );
                    double value = Math.pow(10.0, predicted);
                    if (!(value > 0.0) || !Double.isFinite(value)) {
                        value = fallback;
                    }
                    priorN[i] = Math.max(
                            1.0E-6,
                            Math.min(1.0E6, value)
                    );
                }
            } else if (!scalarPrior.isEmpty()) {
                double value = median(scalarPrior);
                Arrays.fill(priorN, value);
            }
        }

        for (int i = 0; i < n; i++) {
            if (!(priorN[i] > 0.0) || !Double.isFinite(priorN[i])) {
                continue;
            }
            for (int j = 0; j < GRID_LENGTH; j++) {
                logLik[i][j] += priorN[i] * trend[i][j];
            }
        }
    }

    private static double[][] movingAverageTrend(
            double[][] values,
            double[] expression,
            boolean[] interior
    ) {
        int n = values.length;
        double[][] trend = new double[n][GRID_LENGTH];

        List<Integer> retained = new ArrayList<Integer>();
        for (int i = 0; i < n; i++) {
            if (interior[i]) retained.add(i);
        }

        Collections.sort(retained, new Comparator<Integer>() {
            @Override
            public int compare(Integer left, Integer right) {
                int c = Double.compare(expression[left], expression[right]);
                return c != 0 ? c : Integer.compare(left, right);
            }
        });

        int width = Math.max(
                1,
                (int) Math.floor(TREND_SPAN * retained.size())
        );
        width = Math.min(width, retained.size());

        for (int pos = 0; pos < retained.size(); pos++) {
            int start = pos - width / 2;
            if (start < 0) start = 0;
            if (start + width > retained.size()) {
                start = retained.size() - width;
            }
            int end = start + width;

            int target = retained.get(pos);
            for (int col = 0; col < GRID_LENGTH; col++) {
                double sum = 0.0;
                for (int k = start; k < end; k++) {
                    sum += values[retained.get(k)][col];
                }
                trend[target][col] = sum / width;
            }
        }

        for (int i = 0; i < n; i++) {
            if (interior[i]) continue;

            int nearest = retained.get(0);
            double best = Math.abs(expression[i] - expression[nearest]);

            for (int candidate : retained) {
                double distance = Math.abs(
                        expression[i] - expression[candidate]
                );
                if (distance < best) {
                    best = distance;
                    nearest = candidate;
                }
            }
            trend[i] = Arrays.copyOf(trend[nearest], GRID_LENGTH);
        }

        return trend;
    }

    private static double maximizeGrid(double[] y) {
        int k = argMax(y);
        double x = GRID_MIN + k;

        if (k == 0 || k == y.length - 1) {
            return x;
        }

        double left = y[k - 1];
        double center = y[k];
        double right = y[k + 1];
        double denominator = left - 2.0 * center + right;

        if (Math.abs(denominator) < 1.0E-12
                || !Double.isFinite(denominator)) {
            return x;
        }

        double offset = 0.5 * (left - right) / denominator;
        if (!Double.isFinite(offset)) return x;
        offset = Math.max(-1.0, Math.min(1.0, offset));
        return x + offset;
    }

    private static double adjustedProfileLikelihood(
            Event event,
            double precision
    ) {
        if (!(precision > 0.0) || !Double.isFinite(precision)) {
            return Double.NaN;
        }

        int[][] group0 = subset(event.counts, event.groups, 0);
        int[][] group1 = subset(event.counts, event.groups, 1);
        if (group0.length == 0 || group1.length == 0) {
            return Double.NaN;
        }

        DirichletMultinomial.FitResult fit0 =
                DirichletMultinomial.fitOneGroup(group0, precision);
        DirichletMultinomial.FitResult fit1 =
                DirichletMultinomial.fitOneGroup(group1, precision);

        /*
         * Precision profiling uses the best finite optimizer state. Do not
         * discard an event solely because the optimizer stopped on a strict
         * gradient/step criterion.
         */
        if (!Double.isFinite(fit0.logLikelihood)
                || !Double.isFinite(fit1.logLikelihood)) {
            return Double.NaN;
        }

        double adjustment0 = coxReidAdjustment(
                group0,
                fit0.proportions,
                precision
        );
        double adjustment1 = coxReidAdjustment(
                group1,
                fit1.proportions,
                precision
        );

        if (!Double.isFinite(adjustment0)
                || !Double.isFinite(adjustment1)) {
            return Double.NaN;
        }

        return fit0.logLikelihood
                + fit1.logLikelihood
                - adjustment0
                - adjustment1;
    }

    private static double coxReidAdjustment(
            int[][] counts,
            double[] proportions,
            double precision
    ) {
        int features = counts[0].length;
        long[] totals = new long[features];

        for (int[] row : counts) {
            for (int j = 0; j < features; j++) {
                totals[j] += row[j];
            }
        }

        List<Integer> active = new ArrayList<Integer>();
        for (int j = 0; j < features; j++) {
            if (totals[j] > 0L) active.add(j);
        }
        if (active.size() < 2) return Double.NaN;

        int reference = active.get(active.size() - 1);
        double propReference = proportions[reference];
        if (!(propReference > 0.0)
                || !Double.isFinite(propReference)) {
            return Double.NaN;
        }

        int d = active.size() - 1;
        double alphaReference = precision * propReference;

        double common = 0.0;
        for (int[] row : counts) {
            common += trigamma(row[reference] + alphaReference)
                    - trigamma(alphaReference);
        }
        common *= precision * precision;

        double[][] minusH = new double[d][d];
        for (int j = 0; j < d; j++) {
            int feature = active.get(j);
            double p = proportions[feature];
            if (!(p > 0.0) || !Double.isFinite(p)) {
                return Double.NaN;
            }

            double alpha = precision * p;
            double diagonal = 0.0;
            for (int[] row : counts) {
                diagonal += trigamma(row[feature] + alpha)
                        - trigamma(alpha);
            }
            diagonal *= precision * precision;

            for (int k = 0; k < d; k++) {
                double h = common;
                if (j == k) h += diagonal;
                minusH[j][k] = -counts.length * h;
            }
        }

        double logDet = logDeterminantPositiveDefinite(minusH);
        return Double.isFinite(logDet)
                ? 0.5 * logDet
                : Double.NaN;
    }

    private static int[][] subset(
            int[][] counts,
            int[] groups,
            int target
    ) {
        int n = 0;
        for (int group : groups) {
            if (group == target) n++;
        }

        int[][] out = new int[n][];
        int cursor = 0;
        for (int i = 0; i < groups.length; i++) {
            if (groups[i] == target) {
                out[cursor++] = Arrays.copyOf(
                        counts[i],
                        counts[i].length
                );
            }
        }
        return out;
    }

    private static double logDeterminantPositiveDefinite(
            double[][] matrix
    ) {
        int n = matrix.length;
        double[][] l = new double[n][n];
        double logDet = 0.0;

        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                double sum = matrix[i][j];
                for (int k = 0; k < j; k++) {
                    sum -= l[i][k] * l[j][k];
                }

                if (i == j) {
                    if (!(sum > 0.0) || !Double.isFinite(sum)) {
                        return Double.NaN;
                    }
                    l[i][j] = Math.sqrt(sum);
                    logDet += 2.0 * Math.log(l[i][j]);
                } else {
                    l[i][j] = sum / l[j][j];
                }
            }
        }
        return logDet;
    }

    static double trigamma(double x) {
        if (!(x > 0.0) || !Double.isFinite(x)) {
            return Double.NaN;
        }

        double result = 0.0;
        while (x < 8.0) {
            result += 1.0 / (x * x);
            x += 1.0;
        }

        double inv = 1.0 / x;
        double inv2 = inv * inv;
        result += inv
                + 0.5 * inv2
                + inv * inv2 / 6.0
                - inv * inv2 * inv2 / 30.0
                + inv * inv2 * inv2 * inv2 / 42.0
                - inv * inv2 * inv2 * inv2 * inv2 / 30.0
                + 5.0 * inv * Math.pow(inv2, 5) / 66.0;
        return result;
    }

    private static int argMax(double[] values) {
        int best = 0;
        for (int i = 1; i < values.length; i++) {
            if (values[i] > values[best]) best = i;
        }
        return best;
    }

    private static double span(double[] values) {
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        for (double value : values) {
            min = Math.min(min, value);
            max = Math.max(max, value);
        }
        return max - min;
    }

    private static double median(List<Double> values) {
        if (values.isEmpty()) return Double.NaN;
        List<Double> copy = new ArrayList<Double>(values);
        Collections.sort(copy);
        int n = copy.size();
        return n % 2 == 1
                ? copy.get(n / 2)
                : 0.5 * (copy.get(n / 2 - 1)
                + copy.get(n / 2));
    }

    private static double localQuadraticLoess(
            final List<Double> xValues,
            final List<Double> yValues,
            final double x0,
            double span
    ) {
        int n = xValues.size();
        if (n == 0) return Double.NaN;

        Integer[] order = new Integer[n];
        for (int i = 0; i < n; i++) order[i] = i;

        Arrays.sort(order, new Comparator<Integer>() {
            @Override
            public int compare(Integer left, Integer right) {
                return Double.compare(
                        Math.abs(xValues.get(left) - x0),
                        Math.abs(xValues.get(right) - x0)
                );
            }
        });

        int keep = Math.max(3, (int) Math.ceil(span * n));
        keep = Math.min(keep, n);

        double maxDistance = Math.abs(
                xValues.get(order[keep - 1]) - x0
        );
        if (!(maxDistance > 0.0)) {
            return yValues.get(order[0]);
        }

        double s0 = 0.0, s1 = 0.0, s2 = 0.0;
        double s3 = 0.0, s4 = 0.0;
        double t0 = 0.0, t1 = 0.0, t2 = 0.0;

        for (int k = 0; k < keep; k++) {
            int i = order[k];
            double d = xValues.get(i) - x0;
            double u = Math.abs(d) / maxDistance;
            if (u >= 1.0) u = Math.nextDown(1.0);
            double oneMinusCube = 1.0 - u * u * u;
            double w = oneMinusCube
                    * oneMinusCube
                    * oneMinusCube;
            double y = yValues.get(i);
            double d2 = d * d;

            s0 += w;
            s1 += w * d;
            s2 += w * d2;
            s3 += w * d2 * d;
            s4 += w * d2 * d2;

            t0 += w * y;
            t1 += w * d * y;
            t2 += w * d2 * y;
        }

        double[][] a = {
                {s0, s1, s2},
                {s1, s2, s3},
                {s2, s3, s4}
        };
        double[] b = {t0, t1, t2};
        double[] beta = solve3x3(a, b);

        if (beta == null || !Double.isFinite(beta[0])) {
            return s0 > 0.0 ? t0 / s0 : Double.NaN;
        }
        return beta[0];
    }

    private static double[] solve3x3(double[][] a, double[] b) {
        double[][] m = new double[3][4];
        for (int i = 0; i < 3; i++) {
            System.arraycopy(a[i], 0, m[i], 0, 3);
            m[i][3] = b[i];
        }

        for (int col = 0; col < 3; col++) {
            int pivot = col;
            for (int row = col + 1; row < 3; row++) {
                if (Math.abs(m[row][col])
                        > Math.abs(m[pivot][col])) {
                    pivot = row;
                }
            }
            if (Math.abs(m[pivot][col]) < 1.0E-12) {
                return null;
            }
            if (pivot != col) {
                double[] tmp = m[pivot];
                m[pivot] = m[col];
                m[col] = tmp;
            }

            double scale = m[col][col];
            for (int j = col; j < 4; j++) {
                m[col][j] /= scale;
            }

            for (int row = 0; row < 3; row++) {
                if (row == col) continue;
                double factor = m[row][col];
                for (int j = col; j < 4; j++) {
                    m[row][j] -= factor * m[col][j];
                }
            }
        }

        return new double[]{m[0][3], m[1][3], m[2][3]};
    }

    private static long stableHash(String value) {
        long hash = 0xcbf29ce484222325L;
        for (int i = 0; i < value.length(); i++) {
            hash ^= value.charAt(i);
            hash *= 0x100000001b3L;
        }
        return hash;
    }
}
