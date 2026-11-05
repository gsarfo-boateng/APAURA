package org.apaura.analysis;

import org.apaura.model.Strand;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * Condition-blind discovery of coverage-supported candidate PAS breakpoints
 * in a known terminal region.
 *
 * Conventional RNA-seq fragment ends are NOT interpreted as cleavage sites.
 * Candidates are inferred from sustained changes in coverage along transcript
 * orientation and must be reproducible across biological samples.
 */
public final class ShortReadPasDiscovery {

    private ShortReadPasDiscovery() {
    }

    public static final class SampleProfile {
        public final String sample;
        public final String chrom;
        public final Strand strand;
        public final int start0;
        public final int end0;
        public final int[] depth;

        public SampleProfile(
                String sample,
                String chrom,
                Strand strand,
                int start0,
                int end0,
                int[] depth
        ) {
            if (sample == null || sample.isEmpty()) {
                throw new IllegalArgumentException("sample is required");
            }
            if (chrom == null || chrom.isEmpty()) {
                throw new IllegalArgumentException("chrom is required");
            }
            if (strand == null || strand == Strand.UNKNOWN) {
                throw new IllegalArgumentException("known strand is required");
            }
            if (end0 <= start0) {
                throw new IllegalArgumentException("end0 must be > start0");
            }
            if (depth == null || depth.length != end0 - start0) {
                throw new IllegalArgumentException(
                        "depth length must equal end0-start0"
                );
            }

            this.sample = sample;
            this.chrom = chrom;
            this.strand = strand;
            this.start0 = start0;
            this.end0 = end0;
            this.depth = depth.clone();
        }
    }

    public static final class Candidate {
        public final int position0;
        public final int supportingSamples;
        public final double upstreamDepth;
        public final double downstreamDepth;
        public final double coverageDrop;
        public final String status;

        Candidate(
                int position0,
                int supportingSamples,
                double upstreamDepth,
                double downstreamDepth,
                double coverageDrop,
                String status
        ) {
            this.position0 = position0;
            this.supportingSamples = supportingSamples;
            this.upstreamDepth = upstreamDepth;
            this.downstreamDepth = downstreamDepth;
            this.coverageDrop = coverageDrop;
            this.status = status;
        }

        public boolean accepted() {
            return "accepted".equals(status);
        }
    }

    public static final class Result {
        public final List<Candidate> candidates;
        public final List<Candidate> accepted;

        Result(List<Candidate> candidates) {
            List<Candidate> acceptedCandidates =
                    new ArrayList<Candidate>();

            for (Candidate candidate : candidates) {
                if (candidate.accepted()) {
                    acceptedCandidates.add(candidate);
                }
            }

            this.candidates = Collections.unmodifiableList(
                    new ArrayList<Candidate>(candidates)
            );
            this.accepted = Collections.unmodifiableList(
                    acceptedCandidates
            );
        }
    }

    /*
     * Package-visible so the BAM/SAM collector can process one biological
     * sample at a time and retain only compact breakpoint evidence rather
     * than holding genome-wide per-base depth for every sample in memory.
     */
    static final class SampleCall {
        final String sample;
        final int position0;
        final double upstreamDepth;
        final double downstreamDepth;
        final double drop;

        SampleCall(
                String sample,
                int position0,
                double upstreamDepth,
                double downstreamDepth,
                double drop
        ) {
            this.sample = sample;
            this.position0 = position0;
            this.upstreamDepth = upstreamDepth;
            this.downstreamDepth = downstreamDepth;
            this.drop = drop;
        }
    }

    public static Result discover(
            List<SampleProfile> profiles,
            int window,
            double minDepth,
            double minDrop,
            int minReplicates,
            int clusterDistance,
            int minSegmentLength
    ) {
        if (profiles == null || profiles.isEmpty()) {
            return new Result(Collections.<Candidate>emptyList());
        }

        validateParameters(
                window,
                minDepth,
                minDrop,
                minReplicates,
                clusterDistance,
                minSegmentLength
        );
        validateCompatibleProfiles(profiles);

        List<SampleCall> calls = new ArrayList<SampleCall>();

        for (SampleProfile profile : profiles) {
            calls.addAll(callSample(
                    profile,
                    window,
                    minDepth,
                    minDrop,
                    clusterDistance,
                    minSegmentLength
            ));
        }

        return aggregateCalls(
                calls,
                minReplicates,
                clusterDistance
        );
    }

    static List<SampleCall> callSample(
            SampleProfile profile,
            int window,
            double minDepth,
            double minDrop,
            int clusterDistance,
            int minSegmentLength
    ) {
        int[] oriented = orientedDepth(profile);
        int n = oriented.length;

        if (n < 2 * minSegmentLength + 1) {
            return Collections.emptyList();
        }

        double[] prefix = new double[n + 1];
        for (int i = 0; i < n; i++) {
            prefix[i + 1] = prefix[i] + oriented[i];
        }

        List<SampleCall> selected =
                new ArrayList<SampleCall>();

        SampleCall bestInRun = null;
        int previousBoundary = -2;

        for (int boundary = minSegmentLength;
                boundary <= n - minSegmentLength;
                boundary++) {

            int localWindow = Math.min(
                    window,
                    Math.min(boundary, n - boundary)
            );

            if (localWindow < 1) {
                continue;
            }

            double upstream = mean(
                    prefix,
                    boundary - localWindow,
                    boundary
            );

            double downstream = mean(
                    prefix,
                    boundary,
                    boundary + localWindow
            );

            boolean qualifies = false;
            double drop = Double.NaN;

            if (upstream >= minDepth && upstream > 0.0) {
                drop = 1.0 - downstream / upstream;
                qualifies = drop >= minDrop;
            }

            if (!qualifies) {
                if (bestInRun != null) {
                    selected.add(bestInRun);
                    bestInRun = null;
                }
                previousBoundary = -2;
                continue;
            }

            int position0 = profile.strand == Strand.MINUS
                    ? profile.end0 - boundary
                    : profile.start0 + boundary - 1;

            SampleCall candidate = new SampleCall(
                    profile.sample,
                    position0,
                    upstream,
                    downstream,
                    drop
            );

            /*
             * Neighboring qualifying boundaries represent the same broad
             * coverage transition. Collapse each contiguous qualifying run
             * to its strongest breakpoint.
             */
            if (bestInRun == null
                    || boundary != previousBoundary + 1) {

                if (bestInRun != null) {
                    selected.add(bestInRun);
                }

                bestInRun = candidate;

            } else if (betterBreakpoint(candidate, bestInRun)) {
                bestInRun = candidate;
            }

            previousBoundary = boundary;
        }

        if (bestInRun != null) {
            selected.add(bestInRun);
        }

        Collections.sort(
                selected,
                Comparator.comparingInt(call -> call.position0)
        );

        return selected;
    }

    private static boolean betterBreakpoint(
            SampleCall candidate,
            SampleCall current
    ) {
        int comparison = Double.compare(
                candidate.drop,
                current.drop
        );

        if (comparison != 0) {
            return comparison > 0;
        }

        comparison = Double.compare(
                candidate.upstreamDepth,
                current.upstreamDepth
        );

        if (comparison != 0) {
            return comparison > 0;
        }

        /*
         * Deterministic tie-breaker. This does not use read endpoints;
         * it simply resolves equal-scoring coverage boundaries.
         */
        return candidate.position0 < current.position0;
    }

    static Result aggregateCalls(
            List<SampleCall> calls,
            int minReplicates,
            int clusterDistance
    ) {
        if (calls == null || calls.isEmpty()) {
            return new Result(Collections.<Candidate>emptyList());
        }

        List<SampleCall> ordered =
                new ArrayList<SampleCall>(calls);

        Collections.sort(
                ordered,
                Comparator.comparingInt(call -> call.position0)
        );

        List<Candidate> candidates =
                new ArrayList<Candidate>();

        List<SampleCall> cluster =
                new ArrayList<SampleCall>();

        for (SampleCall call : ordered) {
            if (cluster.isEmpty()) {
                cluster.add(call);
                continue;
            }

            int center = medianPosition(cluster);

            if (Math.abs(call.position0 - center)
                    <= clusterDistance) {
                cluster.add(call);
            } else {
                candidates.add(summarize(
                        cluster,
                        minReplicates
                ));
                cluster.clear();
                cluster.add(call);
            }
        }

        if (!cluster.isEmpty()) {
            candidates.add(summarize(
                    cluster,
                    minReplicates
            ));
        }

        return new Result(candidates);
    }

    private static Candidate summarize(
            List<SampleCall> cluster,
            int minReplicates
    ) {
        Set<String> samples = new LinkedHashSet<String>();

        double upstream = 0.0;
        double downstream = 0.0;
        double drop = 0.0;

        for (SampleCall call : cluster) {
            samples.add(call.sample);
            upstream += call.upstreamDepth;
            downstream += call.downstreamDepth;
            drop += call.drop;
        }

        int n = cluster.size();
        int support = samples.size();

        return new Candidate(
                medianPosition(cluster),
                support,
                upstream / n,
                downstream / n,
                drop / n,
                support >= minReplicates
                        ? "accepted"
                        : "insufficient_replicate_support"
        );
    }

    private static int medianPosition(
            List<SampleCall> calls
    ) {
        List<Integer> positions =
                new ArrayList<Integer>();

        for (SampleCall call : calls) {
            positions.add(call.position0);
        }

        Collections.sort(positions);

        return positions.get(
                (positions.size() - 1) / 2
        );
    }

    private static double mean(
            double[] prefix,
            int start,
            int end
    ) {
        return (prefix[end] - prefix[start])
                / (double) (end - start);
    }

    private static int[] orientedDepth(
            SampleProfile profile
    ) {
        if (profile.strand != Strand.MINUS) {
            return profile.depth.clone();
        }

        int[] reversed =
                new int[profile.depth.length];

        for (int i = 0; i < profile.depth.length; i++) {
            reversed[i] =
                    profile.depth[
                            profile.depth.length - 1 - i
                    ];
        }

        return reversed;
    }

    private static void validateParameters(
            int window,
            double minDepth,
            double minDrop,
            int minReplicates,
            int clusterDistance,
            int minSegmentLength
    ) {
        if (window < 1) {
            throw new IllegalArgumentException(
                    "window must be >= 1"
            );
        }

        if (minDepth < 0.0) {
            throw new IllegalArgumentException(
                    "minDepth must be >= 0"
            );
        }

        if (minDrop <= 0.0 || minDrop >= 1.0) {
            throw new IllegalArgumentException(
                    "minDrop must be > 0 and < 1"
            );
        }

        if (minReplicates < 1) {
            throw new IllegalArgumentException(
                    "minReplicates must be >= 1"
            );
        }

        if (clusterDistance < 0) {
            throw new IllegalArgumentException(
                    "clusterDistance must be >= 0"
            );
        }

        if (minSegmentLength < 1) {
            throw new IllegalArgumentException(
                    "minSegmentLength must be >= 1"
            );
        }
    }

    private static void validateCompatibleProfiles(
            List<SampleProfile> profiles
    ) {
        SampleProfile first = profiles.get(0);

        for (SampleProfile profile : profiles) {
            if (!first.chrom.equals(profile.chrom)
                    || first.strand != profile.strand
                    || first.start0 != profile.start0
                    || first.end0 != profile.end0) {

                throw new IllegalArgumentException(
                        "all profiles must describe the same terminal region"
                );
            }
        }
    }
}
