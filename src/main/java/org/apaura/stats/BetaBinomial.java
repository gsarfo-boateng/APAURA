package org.apaura.stats;

/**
 * Replicate-aware beta-binomial likelihood-ratio test for differential site usage.
 * Group labels are encoded as 0 (control), 1 (case), or any other value (ignored).
 */
public final class BetaBinomial {
    private static final double MIN_PROBABILITY = 1.0E-8;
    private static final double MIN_RHO = 1.0E-6;
    private static final double MAX_RHO = 0.95;

    private BetaBinomial() {
    }

    public static final class TestResult {
        public final double pvalue;
        public final double rho;
        public final double p0;
        public final double p1;
        public final double llNull;
        public final double llAlt;

        TestResult(
                double pValue,
                double dispersion,
                double controlProbability,
                double caseProbability,
                double nullLogLikelihood,
                double alternativeLogLikelihood
        ) {
            this.pvalue = pValue;
            this.rho = dispersion;
            this.p0 = controlProbability;
            this.p1 = caseProbability;
            this.llNull = nullLogLikelihood;
            this.llAlt = alternativeLogLikelihood;
        }
    }

    private interface ScalarFunction {
        double evaluate(double x);
    }

    private static final class Fit {
        double logLikelihood;
        double rho;
        double controlProbability;
        double caseProbability;
    }

    private static double logPmf(int successes, int trials, double probability, double rho) {
        double boundedProbability = Math.max(
                MIN_PROBABILITY,
                Math.min(1.0 - MIN_PROBABILITY, probability)
        );
        double boundedRho = Math.max(MIN_RHO, Math.min(MAX_RHO, rho));
        double concentration = 1.0 / boundedRho - 1.0;
        double alpha = boundedProbability * concentration;
        double beta = (1.0 - boundedProbability) * concentration;

        return MathUtil.logChoose(trials, successes)
                + MathUtil.logGamma(successes + alpha)
                + MathUtil.logGamma(trials - successes + beta)
                - MathUtil.logGamma(trials + alpha + beta)
                + MathUtil.logGamma(alpha + beta)
                - MathUtil.logGamma(alpha)
                - MathUtil.logGamma(beta);
    }

    private static double logLikelihood(
            int[] successes,
            int[] totals,
            int[] sampleIndices,
            double probability,
            double rho
    ) {
        double logLikelihood = 0.0;
        for (int sampleIndex : sampleIndices) {
            if (totals[sampleIndex] <= 0) {
                continue;
            }
            logLikelihood += logPmf(
                    successes[sampleIndex],
                    totals[sampleIndex],
                    probability,
                    rho
            );
        }
        return logLikelihood;
    }

    private static double goldenSectionMaximum(
            ScalarFunction function,
            double lower,
            double upper
    ) {
        double goldenRatio = (Math.sqrt(5.0) - 1.0) / 2.0;
        double leftPoint = upper - goldenRatio * (upper - lower);
        double rightPoint = lower + goldenRatio * (upper - lower);
        double leftValue = function.evaluate(leftPoint);
        double rightValue = function.evaluate(rightPoint);

        for (int iteration = 0; iteration < 70; iteration++) {
            if (leftValue > rightValue) {
                upper = rightPoint;
                rightPoint = leftPoint;
                rightValue = leftValue;
                leftPoint = upper - goldenRatio * (upper - lower);
                leftValue = function.evaluate(leftPoint);
            } else {
                lower = leftPoint;
                leftPoint = rightPoint;
                leftValue = rightValue;
                rightPoint = lower + goldenRatio * (upper - lower);
                rightValue = function.evaluate(rightPoint);
            }
        }
        return (lower + upper) / 2.0;
    }

    private static double bestProbability(
            final int[] successes,
            final int[] totals,
            final int[] sampleIndices,
            final double rho
    ) {
        return goldenSectionMaximum(
                new ScalarFunction() {
                    @Override
                    public double evaluate(double probability) {
                        return logLikelihood(
                                successes,
                                totals,
                                sampleIndices,
                                probability,
                                rho
                        );
                    }
                },
                1.0E-6,
                1.0 - 1.0E-6
        );
    }

    private static double rhoFromLogit(double value) {
        return MAX_RHO / (1.0 + Math.exp(-value));
    }

    private static Fit fit(
            final int[] successes,
            final int[] totals,
            final int[] groups,
            final boolean alternative
    ) {
        final int[] allSamples = indices(groups, -1);
        final int[] controlSamples = indices(groups, 0);
        final int[] caseSamples = indices(groups, 1);

        ScalarFunction dispersionObjective = new ScalarFunction() {
            @Override
            public double evaluate(double logitRho) {
                double rho = rhoFromLogit(logitRho);
                if (!alternative) {
                    double probability = bestProbability(
                            successes,
                            totals,
                            allSamples,
                            rho
                    );
                    return logLikelihood(
                            successes,
                            totals,
                            allSamples,
                            probability,
                            rho
                    );
                }

                double controlProbability = bestProbability(
                        successes,
                        totals,
                        controlSamples,
                        rho
                );
                double caseProbability = bestProbability(
                        successes,
                        totals,
                        caseSamples,
                        rho
                );
                return logLikelihood(
                        successes,
                        totals,
                        controlSamples,
                        controlProbability,
                        rho
                ) + logLikelihood(
                        successes,
                        totals,
                        caseSamples,
                        caseProbability,
                        rho
                );
            }
        };

        double bestLogitRho = -8.0;
        double bestLogLikelihood = dispersionObjective.evaluate(bestLogitRho);
        for (int i = 1; i <= 24; i++) {
            double candidateLogitRho = -12.0 + i * (16.0 / 24.0);
            double candidateLogLikelihood = dispersionObjective.evaluate(candidateLogitRho);
            if (candidateLogLikelihood > bestLogLikelihood) {
                bestLogLikelihood = candidateLogLikelihood;
                bestLogitRho = candidateLogitRho;
            }
        }

        double lower = Math.max(-12.0, bestLogitRho - 1.5);
        double upper = Math.min(4.0, bestLogitRho + 1.5);
        double optimizedLogitRho = goldenSectionMaximum(
                dispersionObjective,
                lower,
                upper
        );
        double rho = rhoFromLogit(optimizedLogitRho);

        Fit fit = new Fit();
        fit.rho = rho;
        if (!alternative) {
            fit.controlProbability = bestProbability(successes, totals, allSamples, rho);
            fit.caseProbability = fit.controlProbability;
            fit.logLikelihood = logLikelihood(
                    successes,
                    totals,
                    allSamples,
                    fit.controlProbability,
                    rho
            );
        } else {
            fit.controlProbability = bestProbability(
                    successes,
                    totals,
                    controlSamples,
                    rho
            );
            fit.caseProbability = bestProbability(
                    successes,
                    totals,
                    caseSamples,
                    rho
            );
            fit.logLikelihood = logLikelihood(
                    successes,
                    totals,
                    controlSamples,
                    fit.controlProbability,
                    rho
            ) + logLikelihood(
                    successes,
                    totals,
                    caseSamples,
                    fit.caseProbability,
                    rho
            );
        }
        return fit;
    }

    private static int[] indices(int[] groups, int requestedGroup) {
        int count = 0;
        for (int group : groups) {
            if (matchesGroup(group, requestedGroup)) {
                count++;
            }
        }

        int[] indices = new int[count];
        int next = 0;
        for (int i = 0; i < groups.length; i++) {
            if (matchesGroup(groups[i], requestedGroup)) {
                indices[next++] = i;
            }
        }
        return indices;
    }

    private static boolean matchesGroup(int group, int requestedGroup) {
        if (requestedGroup == -1) {
            return group == 0 || group == 1;
        }
        return group == requestedGroup;
    }

    public static TestResult test(int[] successes, int[] totals, int[] groups) {
        Fit nullFit = fit(successes, totals, groups, false);
        Fit alternativeFit = fit(successes, totals, groups, true);

        double likelihoodRatio = Math.max(
                0.0,
                2.0 * (alternativeFit.logLikelihood - nullFit.logLikelihood)
        );
        double pValue = MathUtil.chiSquare1Survival(likelihoodRatio);

        return new TestResult(
                pValue,
                alternativeFit.rho,
                alternativeFit.controlProbability,
                alternativeFit.caseProbability,
                nullFit.logLikelihood,
                alternativeFit.logLikelihood
        );
    }
}
