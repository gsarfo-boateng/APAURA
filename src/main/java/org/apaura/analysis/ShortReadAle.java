package org.apaura.analysis;

import org.apaura.model.CoverageSegment;
import org.apaura.model.Interval;
import org.apaura.model.Sample;
import org.apaura.model.Strand;
import org.apaura.model.TranscriptModel;
import org.apaura.stats.BetaBinomial;
import org.apaura.stats.FisherExact;
import org.apaura.stats.MultipleTesting;
import org.apaura.stats.ModeratedDirichletMultinomial;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/**
 * Annotation-guided alternative-last-exon (ALE) quantification for conventional
 * short-read RNA-seq.
 *
 * <p>Terminal-exon branches are obtained from the same terminal-family
 * reconciliation used by ArchitectureCatalog and ShortReadGroupBuilder.
 * Evidence is branch-specific: an alignment supports a branch when it overlaps
 * branch-common sequence unique to that terminal branch or contains the splice boundary that
 * enters that branch. Alignments supporting more than one branch of the same
 * gene are treated as ambiguous and are not counted for ALE usage.</p>
 */
public final class ShortReadAle {

    private static final int DEFAULT_MIN_INFORMATIVE_REPLICATES = 3;

    private ShortReadAle() {
        // Utility class.
    }

    public static final class Branch {
        public final String id;
        public final String geneId;
        public final String geneName;
        public final String biotype;
        public final String chrom;
        public final Strand strand;
        public final int rank;
        public final int anchor0;
        public final int representativePas0;
        public final List<Integer> annotatedPas0;
        public final List<CoverageSegment> uniqueRegions;
        public final int availableUniqueBases;
        public final int evidenceBases;
        public final double normalizedPosition;

        Branch(
                String id,
                String geneId,
                String geneName,
                String biotype,
                String chrom,
                Strand strand,
                int rank,
                int anchor0,
                int representativePas0,
                List<Integer> annotatedPas0,
                List<CoverageSegment> uniqueRegions,
                int availableUniqueBases,
                int evidenceBases,
                double normalizedPosition
        ) {
            this.id = id;
            this.geneId = geneId;
            this.geneName = geneName;
            this.biotype = biotype;
            this.chrom = chrom;
            this.strand = strand;
            this.rank = rank;
            this.anchor0 = anchor0;
            this.representativePas0 = representativePas0;
            this.annotatedPas0 = Collections.unmodifiableList(
                    new ArrayList<Integer>(annotatedPas0)
            );
            this.uniqueRegions = Collections.unmodifiableList(
                    new ArrayList<CoverageSegment>(uniqueRegions)
            );
            this.availableUniqueBases = availableUniqueBases;
            this.evidenceBases = evidenceBases;
            this.normalizedPosition = normalizedPosition;
        }
    }

    public static final class Group {
        public final String id;
        public final String geneId;
        public final String geneName;
        public final String biotype;
        public final String chrom;
        public final Strand strand;
        public final List<Branch> branches;

        Group(
                String id,
                String geneId,
                String geneName,
                String biotype,
                String chrom,
                Strand strand,
                List<Branch> branches
        ) {
            this.id = id;
            this.geneId = geneId;
            this.geneName = geneName;
            this.biotype = biotype;
            this.chrom = chrom;
            this.strand = strand;
            this.branches = Collections.unmodifiableList(
                    new ArrayList<Branch>(branches)
            );
        }
    }

    public static final class Quantification {
        public final Group group;
        public final int[] branchReads;
        public final double[] branchUsage;
        public final int totalEvidenceReads;

        Quantification(Group group, int[] branchReads) {
            this.group = group;
            this.branchReads = Arrays.copyOf(branchReads, branchReads.length);
            int total = 0;
            for (int count : branchReads) {
                total += count;
            }
            this.totalEvidenceReads = total;
            this.branchUsage = new double[branchReads.length];
            if (total > 0) {
                for (int index = 0; index < branchReads.length; index++) {
                    this.branchUsage[index] = branchReads[index] / (double) total;
                }
            }
        }
    }

    public static final class BranchResult {
        public Group group;
        public Branch branch;
        public double controlUsage = Double.NaN;
        public double caseUsage = Double.NaN;
        public double delta = Double.NaN;
        public double controlUsageSd = Double.NaN;
        public double caseUsageSd = Double.NaN;
        public int controlTotalReplicates;
        public int caseTotalReplicates;
        public int controlAnyEvidence;
        public int caseAnyEvidence;
        public int controlInformativeReplicates;
        public int caseInformativeReplicates;
        public int controlEvidence;
        public int caseEvidence;
        public double p = Double.NaN;
        public double q = Double.NaN;
        public double rho = Double.NaN;
        public String method = "not_tested";
        public String status = "insufficient_evidence";
    }

    public static final class GeneResult {
        public Group group;
        public double controlWpi = Double.NaN;
        public double caseWpi = Double.NaN;
        public double deltaWpi = Double.NaN;
        public double controlWpiSd = Double.NaN;
        public double caseWpiSd = Double.NaN;
        public double switchIndex = Double.NaN;
        public double minPvalue = Double.NaN;
        public double minFdr = Double.NaN;
        public double eventPvalue = Double.NaN;
        public double eventFdr = Double.NaN;
        public double eventPrecision = Double.NaN;
        public String eventMethod = "not_tested";
        public double minBranchPvalue = Double.NaN;
        public double minBranchFdr = Double.NaN;
        public int controlTotalReplicates;
        public int caseTotalReplicates;
        public int controlAnyEvidence;
        public int caseAnyEvidence;
        public int controlInformativeReplicates;
        public int caseInformativeReplicates;
        public boolean aleChange;
        public String observedDirection = "not_applicable";
        public String aleDirection = "not_applicable";
        public String status = "insufficient_evidence";
    }

    public static final class Results {
        public final List<BranchResult> branches = new ArrayList<BranchResult>();
        public final List<GeneResult> genes = new ArrayList<GeneResult>();
    }

    public static List<Group> build(
            List<TranscriptModel> transcripts,
            int anchorTolerance,
            int pasMergeDistance,
            int minUniqueBases
    ) {
        Map<String, List<TranscriptModel>> byGene = groupByGene(transcripts);
        List<Group> groups = new ArrayList<Group>();

        for (List<TranscriptModel> geneTranscripts : byGene.values()) {
            List<TranscriptModel> terminalTranscripts =
                    ShortReadIpa.terminalArchitectureTranscripts(geneTranscripts);
            if (terminalTranscripts.isEmpty()) {
                continue;
            }
            TerminalFamilyReconciler.Result reconciliation =
                    TerminalFamilyReconciler.reconcile(
                            terminalTranscripts,
                            anchorTolerance,
                            pasMergeDistance
                    );
            if (reconciliation.families.size() < 2) {
                continue;
            }
            groups.add(toGroup(reconciliation.families, minUniqueBases));
        }

        Collections.sort(groups, new Comparator<Group>() {
            @Override
            public int compare(Group left, Group right) {
                int c = left.chrom.compareTo(right.chrom);
                return c != 0 ? c : left.geneId.compareTo(right.geneId);
            }
        });
        return Collections.unmodifiableList(groups);
    }

    private static Group toGroup(
            List<TerminalFamilyReconciler.Family> families,
            int minUniqueBases
    ) {
        TerminalFamilyReconciler.Family first = families.get(0);
        List<FamilyView> views = new ArrayList<FamilyView>();

        for (TerminalFamilyReconciler.Family family : families) {
            /*
             * ALE evidence must be independent of tandem-PAS choice within a
             * branch. Therefore use only the branch-common terminal segment:
             * anchor -> most proximal PAS. Distal TUTR-only sequence is not
             * allowed to contribute to branch usage.
             */
            int proximalPas0 = family.pasPositions.get(0);
            int start0;
            int end0;
            if (family.strand == Strand.MINUS) {
                start0 = proximalPas0;
                end0 = family.anchor0;
            } else {
                start0 = family.anchor0;
                end0 = proximalPas0 + 1;
            }
            if (end0 <= start0) {
                start0 = Math.min(family.anchor0, proximalPas0);
                end0 = Math.max(family.anchor0, proximalPas0 + 1);
            }

            FamilyView view = new FamilyView();
            view.family = family;
            view.span = new Interval(start0, end0);
            view.representativePas0 = median(family.pasPositions);
            views.add(view);
        }

        Collections.sort(views, new Comparator<FamilyView>() {
            @Override
            public int compare(FamilyView left, FamilyView right) {
                int c = Integer.compare(left.representativePas0, right.representativePas0);
                return first.strand == Strand.MINUS ? -c : c;
            }
        });

        /*
         * Raw branch-unique sequence can differ greatly in length between ALE
         * branches. Counting reads across unequal regions would bias branch
         * usage toward the branch with the longer unique footprint. Build the
         * unique pieces first, then use the SAME number of unique bases for
         * every branch in the gene. If any branch lacks the minimum amount of
         * unique sequence, disable exonic coverage for the whole ALE group and
         * rely on branch-specific splice-junction evidence only.
         */
        int minimumAvailableUniqueBases = Integer.MAX_VALUE;
        for (int index = 0; index < views.size(); index++) {
            FamilyView view = views.get(index);
            List<Interval> otherSpans = new ArrayList<Interval>();
            for (int other = 0; other < views.size(); other++) {
                if (other != index) otherSpans.add(views.get(other).span);
            }
            view.unique = subtract(view.span, otherSpans);
            view.availableUniqueBases = totalLength(view.unique);
            minimumAvailableUniqueBases = Math.min(
                    minimumAvailableUniqueBases,
                    view.availableUniqueBases
            );
        }

        int evidenceBases = minimumAvailableUniqueBases >= minUniqueBases
                ? minimumAvailableUniqueBases
                : 0;

        List<Branch> branches = new ArrayList<Branch>();
        int n = views.size();
        for (int index = 0; index < n; index++) {
            FamilyView view = views.get(index);
            List<CoverageSegment> regions = selectEvidenceRegions(
                    view.unique,
                    view.family.anchor0,
                    view.family.strand,
                    evidenceBases
            );

            double normalizedPosition = n == 1
                    ? 0.5
                    : index / (double) (n - 1);
            String id = first.geneId
                    + "|ALEB"
                    + String.format(Locale.ROOT, "%03d", index + 1);

            branches.add(new Branch(
                    id,
                    view.family.geneId,
                    view.family.geneName,
                    view.family.biotype,
                    view.family.chrom,
                    view.family.strand,
                    index + 1,
                    view.family.anchor0,
                    view.representativePas0,
                    view.family.pasPositions,
                    regions,
                    view.availableUniqueBases,
                    evidenceBases,
                    normalizedPosition
            ));
        }

        return new Group(
                first.geneId + "|ALE",
                first.geneId,
                first.geneName,
                preferredBiotype(families),
                first.chrom,
                first.strand,
                branches
        );
    }

    private static final class FamilyView {
        TerminalFamilyReconciler.Family family;
        Interval span;
        int representativePas0;
        List<Interval> unique;
        int availableUniqueBases;
    }

    private static int totalLength(List<Interval> intervals) {
        int total = 0;
        for (Interval interval : intervals) {
            total += Math.max(0, interval.end0 - interval.start0);
        }
        return total;
    }

    private static List<CoverageSegment> selectEvidenceRegions(
            List<Interval> unique,
            int anchor0,
            Strand strand,
            int targetBases
    ) {
        if (targetBases <= 0 || unique.isEmpty()) {
            return Collections.emptyList();
        }

        List<Interval> ordered = new ArrayList<Interval>(unique);
        Collections.sort(ordered, new Comparator<Interval>() {
            @Override
            public int compare(Interval left, Interval right) {
                int leftDistance = distanceToAnchor(left, anchor0);
                int rightDistance = distanceToAnchor(right, anchor0);
                int c = Integer.compare(leftDistance, rightDistance);
                if (c != 0) return c;
                return strand == Strand.MINUS
                        ? Integer.compare(right.end0, left.end0)
                        : Integer.compare(left.start0, right.start0);
            }
        });

        List<CoverageSegment> regions = new ArrayList<CoverageSegment>();
        int remaining = targetBases;
        int regionNumber = 0;
        for (Interval interval : ordered) {
            if (remaining <= 0) break;
            int length = interval.end0 - interval.start0;
            if (length <= 0) continue;
            int take = Math.min(length, remaining);
            int start0;
            int end0;
            if (strand == Strand.MINUS) {
                end0 = interval.end0;
                start0 = end0 - take;
            } else {
                start0 = interval.start0;
                end0 = start0 + take;
            }
            regions.add(new CoverageSegment(
                    regionNumber++, start0, end0, "ALE_UNIQUE"
            ));
            remaining -= take;
        }
        if (remaining != 0) {
            throw new IllegalStateException(
                    "Internal ALE evidence error: requested "
                            + targetBases
                            + " unique bases but only selected "
                            + (targetBases - remaining)
            );
        }
        return regions;
    }

    private static int distanceToAnchor(Interval interval, int anchor0) {
        if (anchor0 < interval.start0) return interval.start0 - anchor0;
        if (anchor0 > interval.end0) return anchor0 - interval.end0;
        return 0;
    }

    private static List<Interval> subtract(Interval source, List<Interval> masks) {
        List<Interval> pieces = new ArrayList<Interval>();
        pieces.add(source);
        for (Interval mask : masks) {
            List<Interval> next = new ArrayList<Interval>();
            for (Interval piece : pieces) {
                if (mask.end0 <= piece.start0 || mask.start0 >= piece.end0) {
                    next.add(piece);
                    continue;
                }
                if (mask.start0 > piece.start0) {
                    next.add(new Interval(piece.start0, Math.min(mask.start0, piece.end0)));
                }
                if (mask.end0 < piece.end0) {
                    next.add(new Interval(Math.max(mask.end0, piece.start0), piece.end0));
                }
            }
            pieces = next;
            if (pieces.isEmpty()) {
                break;
            }
        }
        return pieces;
    }

    static Quantification quantify(Group group, int[] branchReads) {
        return new Quantification(group, branchReads);
    }

    public static Results differential(
            List<Sample> samples,
            Map<String, List<Quantification>> bySample,
            String control,
            String caseCondition,
            int minSampleEvidence,
            int minInformativeReplicates,
            double minBranchDelta,
            double minSwitch,
            double fdr,
            double wpiThreshold
    ) {
        if (minInformativeReplicates < 1) {
            throw new IllegalArgumentException("minInformativeReplicates must be >= 1");
        }
        Results results = new Results();
        if (samples.isEmpty()) {
            return results;
        }
        List<Quantification> first = bySample.get(samples.get(0).name);
        if (first == null) {
            return results;
        }

        Map<String, List<BranchResult>> byGroup =
                new LinkedHashMap<String, List<BranchResult>>();

        for (int groupIndex = 0; groupIndex < first.size(); groupIndex++) {
            Group group = first.get(groupIndex).group;
            for (int branchIndex = 0; branchIndex < group.branches.size(); branchIndex++) {
                BranchResult result = testBranch(
                        samples,
                        bySample,
                        groupIndex,
                        branchIndex,
                        control,
                        caseCondition,
                        minSampleEvidence,
                        minInformativeReplicates
                );
                results.branches.add(result);
                byGroup.computeIfAbsent(
                        group.id,
                        key -> new ArrayList<BranchResult>()
                ).add(result);
            }
        }

        applyFdr(results.branches);

        for (int groupIndex = 0; groupIndex < first.size(); groupIndex++) {
            results.genes.add(summarizeGene(
                    samples,
                    bySample,
                    groupIndex,
                    byGroup.get(first.get(groupIndex).group.id),
                    control,
                    caseCondition,
                    minSampleEvidence,
                    minInformativeReplicates,
                    minBranchDelta,
                    minSwitch,
                    fdr,
                    wpiThreshold
            ));
        }

        applyModeratedDmToGenes(
                samples,
                bySample,
                control,
                caseCondition,
                minSampleEvidence,
                results.genes
        );
        applyGeneFdr(results.genes);
        classifyGenes(results.genes, minSwitch, fdr, wpiThreshold);

        Collections.sort(results.branches, new Comparator<BranchResult>() {
            @Override
            public int compare(BranchResult left, BranchResult right) {
                return compareFinite(left.q, right.q);
            }
        });
        Collections.sort(results.genes, new Comparator<GeneResult>() {
            @Override
            public int compare(GeneResult left, GeneResult right) {
                return compareFinite(left.minFdr, right.minFdr);
            }
        });
        return results;
    }

    public static Results differential(
            List<Sample> samples,
            Map<String, List<Quantification>> bySample,
            String control,
            String caseCondition,
            int minSampleEvidence,
            double minBranchDelta,
            double minSwitch,
            double fdr,
            double wpiThreshold
    ) {
        return differential(
                samples,
                bySample,
                control,
                caseCondition,
                minSampleEvidence,
                DEFAULT_MIN_INFORMATIVE_REPLICATES,
                minBranchDelta,
                minSwitch,
                fdr,
                wpiThreshold
        );
    }

    private static BranchResult testBranch(
            List<Sample> samples,
            Map<String, List<Quantification>> bySample,
            int groupIndex,
            int branchIndex,
            String control,
            String caseCondition,
            int minSampleEvidence,
            int minInformativeReplicates
    ) {
        Group group = bySample.get(samples.get(0).name).get(groupIndex).group;
        Branch branch = group.branches.get(branchIndex);
        BranchResult result = new BranchResult();
        result.group = group;
        result.branch = branch;

        int[] branchCounts = new int[samples.size()];
        int[] totals = new int[samples.size()];
        int[] labels = new int[samples.size()];
        Arrays.fill(labels, -2);
        List<Double> controlUsage = new ArrayList<Double>();
        List<Double> caseUsage = new ArrayList<Double>();

        for (int index = 0; index < samples.size(); index++) {
            Sample sample = samples.get(index);
            Quantification q = bySample.get(sample.name).get(groupIndex);
            totals[index] = q.totalEvidenceReads;
            branchCounts[index] = q.branchReads[branchIndex];

            boolean isControl = sample.condition.equals(control);
            boolean isCase = sample.condition.equals(caseCondition);
            if (!isControl && !isCase) {
                continue;
            }

            if (isControl) {
                result.controlTotalReplicates++;
                if (q.totalEvidenceReads > 0) result.controlAnyEvidence++;
            } else {
                result.caseTotalReplicates++;
                if (q.totalEvidenceReads > 0) result.caseAnyEvidence++;
            }

            if (q.totalEvidenceReads < minSampleEvidence) {
                continue;
            }

            if (isControl) {
                labels[index] = 0;
                result.controlInformativeReplicates++;
                result.controlEvidence += q.totalEvidenceReads;
                controlUsage.add(q.branchUsage[branchIndex]);
            } else {
                labels[index] = 1;
                result.caseInformativeReplicates++;
                result.caseEvidence += q.totalEvidenceReads;
                caseUsage.add(q.branchUsage[branchIndex]);
            }
        }

        result.controlUsage = mean(controlUsage);
        result.caseUsage = mean(caseUsage);
        result.controlUsageSd = sampleSd(controlUsage);
        result.caseUsageSd = sampleSd(caseUsage);
        result.delta = finiteDifference(result.caseUsage, result.controlUsage);

        if (result.controlInformativeReplicates == 0
                || result.caseInformativeReplicates == 0) {
            result.status = "insufficient_evidence";
            return result;
        }

        boolean unreplicated = result.controlTotalReplicates == 1
                && result.caseTotalReplicates == 1;
        if (unreplicated) {
            int[] table = pooledTwoByTwo(branchCounts, totals, labels);
            result.p = FisherExact.twoSided(table[0], table[1], table[2], table[3]);
            result.method = "ale-Fisher-exact-unreplicated";
            result.status = "ok";
            return result;
        }

        if (result.controlInformativeReplicates < minInformativeReplicates
                || result.caseInformativeReplicates < minInformativeReplicates) {
            result.status = "insufficient_replicate_support";
            return result;
        }

        BetaBinomial.TestResult test = BetaBinomial.test(branchCounts, totals, labels);
        result.p = test.pvalue;
        result.rho = test.rho;
        result.method = "ale-beta-binomial-LRT";
        result.status = "ok";
        return result;
    }

    private static GeneResult summarizeGene(
            List<Sample> samples,
            Map<String, List<Quantification>> bySample,
            int groupIndex,
            List<BranchResult> branchResults,
            String control,
            String caseCondition,
            int minSampleEvidence,
            int minInformativeReplicates,
            double minBranchDelta,
            double minSwitch,
            double fdr,
            double wpiThreshold
    ) {
        Group group = bySample.get(samples.get(0).name).get(groupIndex).group;
        GeneResult result = new GeneResult();
        result.group = group;
        double[] controlMean = new double[group.branches.size()];
        double[] caseMean = new double[group.branches.size()];
        List<Double> controlWpi = new ArrayList<Double>();
        List<Double> caseWpi = new ArrayList<Double>();

        for (Sample sample : samples) {
            Quantification q = bySample.get(sample.name).get(groupIndex);
            boolean isControl = sample.condition.equals(control);
            boolean isCase = sample.condition.equals(caseCondition);
            if (!isControl && !isCase) continue;

            if (isControl) {
                result.controlTotalReplicates++;
                if (q.totalEvidenceReads > 0) result.controlAnyEvidence++;
            } else {
                result.caseTotalReplicates++;
                if (q.totalEvidenceReads > 0) result.caseAnyEvidence++;
            }

            if (q.totalEvidenceReads < minSampleEvidence) continue;

            double wpi = weightedPositionIndex(group, q.branchUsage);
            if (isControl) {
                result.controlInformativeReplicates++;
                controlWpi.add(wpi);
                addInPlace(controlMean, q.branchUsage);
            } else {
                result.caseInformativeReplicates++;
                caseWpi.add(wpi);
                addInPlace(caseMean, q.branchUsage);
            }
        }

        if (!controlWpi.isEmpty()) {
            result.controlWpi = mean(controlWpi);
            result.controlWpiSd = sampleSd(controlWpi);
            averageInPlace(controlMean, controlWpi.size());
        }
        if (!caseWpi.isEmpty()) {
            result.caseWpi = mean(caseWpi);
            result.caseWpiSd = sampleSd(caseWpi);
            averageInPlace(caseMean, caseWpi.size());
        }
        result.deltaWpi = finiteDifference(result.caseWpi, result.controlWpi);
        if (!controlWpi.isEmpty() && !caseWpi.isEmpty()) {
            result.switchIndex = totalVariation(controlMean, caseMean);
            result.observedDirection = result.switchIndex >= minSwitch
                    ? observedDirection(result.deltaWpi, wpiThreshold)
                    : "stable";
        }

        if (result.controlInformativeReplicates == 0
                || result.caseInformativeReplicates == 0) {
            result.status = "insufficient_evidence";
            return result;
        }
        boolean unreplicated = result.controlTotalReplicates == 1
                && result.caseTotalReplicates == 1;
        if (!unreplicated
                && (result.controlInformativeReplicates < minInformativeReplicates
                || result.caseInformativeReplicates < minInformativeReplicates)) {
            result.status = "insufficient_replicate_support";
            return result;
        }

        result.status = "ok";
        int testedBranches = 0;
        if (branchResults != null) {
            for (BranchResult branch : branchResults) {
                if ("ok".equals(branch.status) && Double.isFinite(branch.p)) {
                    testedBranches++;
                    if (!Double.isFinite(result.minBranchPvalue)
                            || branch.p < result.minBranchPvalue) {
                        result.minBranchPvalue = branch.p;
                    }
                }
                if (Double.isFinite(branch.q)
                        && (!Double.isFinite(result.minBranchFdr)
                        || branch.q < result.minBranchFdr)) {
                    result.minBranchFdr = branch.q;
                }
            }
        }

        if (unreplicated) {
            if (testedBranches == 0 || !Double.isFinite(result.minBranchPvalue)) {
                result.status = "insufficient_evidence";
                return result;
            }
            result.eventPvalue = Math.min(
                    1.0,
                    result.minBranchPvalue * testedBranches
            );
            result.eventMethod = "bonferroni-min-branch-unreplicated";
            result.minPvalue = result.eventPvalue;
        }

        return result;
    }

    private static void applyModeratedDmToGenes(
            List<Sample> samples,
            Map<String, List<Quantification>> bySample,
            String control,
            String caseCondition,
            int minSampleEvidence,
            List<GeneResult> genes
    ) {
        List<ModeratedDirichletMultinomial.Event> events =
                new ArrayList<ModeratedDirichletMultinomial.Event>();
        Map<String, GeneResult> byId =
                new LinkedHashMap<String, GeneResult>();

        for (int groupIndex = 0; groupIndex < genes.size(); groupIndex++) {
            GeneResult result = genes.get(groupIndex);
            if (!"ok".equals(result.status)) continue;

            boolean unreplicated = result.controlTotalReplicates == 1
                    && result.caseTotalReplicates == 1;
            if (unreplicated) continue;

            List<int[]> rows = new ArrayList<int[]>();
            List<Integer> labels = new ArrayList<Integer>();

            for (Sample sample : samples) {
                boolean isControl = sample.condition.equals(control);
                boolean isCase = sample.condition.equals(caseCondition);
                if (!isControl && !isCase) continue;

                Quantification q = bySample.get(sample.name).get(groupIndex);
                if (q.totalEvidenceReads < minSampleEvidence) continue;

                rows.add(Arrays.copyOf(q.branchReads, q.branchReads.length));
                labels.add(isControl ? 0 : 1);
            }

            if (rows.isEmpty()) continue;

            int[][] counts = new int[rows.size()][];
            int[] groupLabels = new int[rows.size()];
            for (int i = 0; i < rows.size(); i++) {
                counts[i] = rows.get(i);
                groupLabels[i] = labels.get(i);
            }

            result.eventPvalue = Double.NaN;
            result.eventFdr = Double.NaN;
            result.minPvalue = Double.NaN;
            result.minFdr = Double.NaN;
            result.eventMethod = "moderated-DM-LRT";

            events.add(new ModeratedDirichletMultinomial.Event(
                    result.group.id,
                    counts,
                    groupLabels
            ));
            byId.put(result.group.id, result);
        }

        ModeratedDirichletMultinomial.BatchResult batch =
                ModeratedDirichletMultinomial.test(events);

        for (ModeratedDirichletMultinomial.Result dm : batch.results) {
            GeneResult result = byId.get(dm.id);
            if (result == null) continue;

            if (dm.converged && Double.isFinite(dm.pvalue)) {
                result.eventPvalue = dm.pvalue;
                result.minPvalue = dm.pvalue;
                result.eventPrecision = dm.precision;
            } else {
                result.status = "dm_fit_failed";
                result.eventMethod = "moderated-DM-LRT-failed";
            }
        }
    }

    private static void applyGeneFdr(List<GeneResult> genes) {
        List<Integer> indices = new ArrayList<Integer>();
        List<Double> pvalues = new ArrayList<Double>();

        for (int i = 0; i < genes.size(); i++) {
            if (Double.isFinite(genes.get(i).eventPvalue)) {
                indices.add(i);
                pvalues.add(genes.get(i).eventPvalue);
            }
        }
        if (indices.isEmpty()) return;

        double[] p = new double[pvalues.size()];
        for (int i = 0; i < p.length; i++) p[i] = pvalues.get(i);
        double[] q = MultipleTesting.bh(p);

        for (int i = 0; i < q.length; i++) {
            GeneResult result = genes.get(indices.get(i));
            result.eventFdr = q[i];
            result.minFdr = q[i];
        }
    }

    private static void classifyGenes(
            List<GeneResult> genes,
            double minSwitch,
            double fdr,
            double wpiThreshold
    ) {
        for (GeneResult result : genes) {
            if (!"ok".equals(result.status)) continue;

            result.aleChange =
                    Double.isFinite(result.eventFdr)
                    && result.eventFdr <= fdr
                    && Double.isFinite(result.switchIndex)
                    && result.switchIndex >= minSwitch;

            result.aleDirection = result.aleChange
                    ? inferredDirection(result.deltaWpi, wpiThreshold)
                    : "stable";
        }
    }

    private static double weightedPositionIndex(Group group, double[] usage) {
        double wpi = 0.0;
        for (int index = 0; index < group.branches.size(); index++) {
            wpi += usage[index] * group.branches.get(index).normalizedPosition;
        }
        return wpi;
    }

    private static String observedDirection(double deltaWpi, double threshold) {
        if (!Double.isFinite(deltaWpi)) return "not_applicable";
        if (deltaWpi >= threshold) return "distal_switch";
        if (deltaWpi <= -threshold) return "proximal_switch";
        return "switch";
    }

    private static String inferredDirection(double deltaWpi, double threshold) {
        return observedDirection(deltaWpi, threshold);
    }

    private static int[] pooledTwoByTwo(int[] counts, int[] totals, int[] labels) {
        int controlSuccess = 0;
        int controlOther = 0;
        int caseSuccess = 0;
        int caseOther = 0;
        for (int index = 0; index < labels.length; index++) {
            if (labels[index] == 0) {
                controlSuccess += counts[index];
                controlOther += totals[index] - counts[index];
            } else if (labels[index] == 1) {
                caseSuccess += counts[index];
                caseOther += totals[index] - counts[index];
            }
        }
        return new int[]{controlSuccess, controlOther, caseSuccess, caseOther};
    }

    private static void applyFdr(List<BranchResult> results) {
        List<Integer> indices = new ArrayList<Integer>();
        List<Double> pvalues = new ArrayList<Double>();
        for (int index = 0; index < results.size(); index++) {
            if (Double.isFinite(results.get(index).p)) {
                indices.add(index);
                pvalues.add(results.get(index).p);
            }
        }
        if (indices.isEmpty()) return;
        double[] p = new double[pvalues.size()];
        for (int index = 0; index < p.length; index++) p[index] = pvalues.get(index);
        double[] q = MultipleTesting.bh(p);
        for (int index = 0; index < q.length; index++) {
            results.get(indices.get(index)).q = q[index];
        }
    }

    public static void writeAnnotation(File file, List<Group> groups) throws IOException {
        try (PrintWriter writer = writer(file)) {
            writer.println(
                    "gene_id\tgene_name\tbiotype\tchrom\tstrand\tbranch_id\tbranch_rank\t"
                            + "anchor0\trepresentative_pas_1based\tn_annotated_pas\t"
                            + "annotated_pas_positions_1based\tavailable_unique_bases\t"
                            + "evidence_bases\tn_unique_regions"
            );
            for (Group group : groups) {
                for (Branch branch : group.branches) {
                    writer.println(
                            group.geneId + "\t"
                                    + group.geneName + "\t"
                                    + group.biotype + "\t"
                                    + group.chrom + "\t"
                                    + group.strand.symbol() + "\t"
                                    + branch.id + "\t"
                                    + branch.rank + "\t"
                                    + branch.anchor0 + "\t"
                                    + (branch.representativePas0 + 1) + "\t"
                                    + branch.annotatedPas0.size() + "\t"
                                    + positions1(branch.annotatedPas0) + "\t"
                                    + branch.availableUniqueBases + "\t"
                                    + branch.evidenceBases + "\t"
                                    + branch.uniqueRegions.size()
                    );
                }
            }
        }
    }

    public static void writeUsage(
            File file,
            List<Sample> samples,
            Map<String, List<Quantification>> bySample
    ) throws IOException {
        try (PrintWriter writer = writer(file)) {
            writer.println(
                    "sample\tcondition\tgene_id\tgene_name\tbranch_id\tbranch_rank\t"
                            + "branch_reads\ttotal_branch_evidence\tbranch_usage"
            );
            for (Sample sample : samples) {
                List<Quantification> qs = bySample.get(sample.name);
                if (qs == null) continue;
                for (Quantification q : qs) {
                    for (int index = 0; index < q.group.branches.size(); index++) {
                        Branch branch = q.group.branches.get(index);
                        writer.println(
                                sample.name + "\t"
                                        + sample.condition + "\t"
                                        + q.group.geneId + "\t"
                                        + q.group.geneName + "\t"
                                        + branch.id + "\t"
                                        + branch.rank + "\t"
                                        + q.branchReads[index] + "\t"
                                        + q.totalEvidenceReads + "\t"
                                        + fmt(q.branchUsage[index], "%.8f")
                        );
                    }
                }
            }
        }
    }

    public static void writeBranches(File file, List<BranchResult> results) throws IOException {
        try (PrintWriter writer = writer(file)) {
            writer.println(
                    "gene_id\tgene_name\tbranch_id\tbranch_rank\tchrom\tstrand\t"
                            + "representative_pas_1based\tcontrol_usage\tcase_usage\tdelta_usage\t"
                            + "control_usage_sd\tcase_usage_sd\tcontrol_n_total\tcase_n_total\t"
                            + "control_n_any_evidence\tcase_n_any_evidence\tcontrol_n_informative\t"
                            + "case_n_informative\tcontrol_evidence\tcase_evidence\tpvalue\tfdr\t"
                            + "dispersion_rho\tmethod\tstatus"
            );
            for (BranchResult result : results) {
                writer.println(
                        result.group.geneId + "\t"
                                + result.group.geneName + "\t"
                                + result.branch.id + "\t"
                                + result.branch.rank + "\t"
                                + result.group.chrom + "\t"
                                + result.group.strand.symbol() + "\t"
                                + (result.branch.representativePas0 + 1) + "\t"
                                + fmt(result.controlUsage, "%.8f") + "\t"
                                + fmt(result.caseUsage, "%.8f") + "\t"
                                + fmt(result.delta, "%.8f") + "\t"
                                + fmt(result.controlUsageSd, "%.8f") + "\t"
                                + fmt(result.caseUsageSd, "%.8f") + "\t"
                                + result.controlTotalReplicates + "\t"
                                + result.caseTotalReplicates + "\t"
                                + result.controlAnyEvidence + "\t"
                                + result.caseAnyEvidence + "\t"
                                + result.controlInformativeReplicates + "\t"
                                + result.caseInformativeReplicates + "\t"
                                + result.controlEvidence + "\t"
                                + result.caseEvidence + "\t"
                                + fmt(result.p, "%.10g") + "\t"
                                + fmt(result.q, "%.10g") + "\t"
                                + fmt(result.rho, "%.8f") + "\t"
                                + result.method + "\t"
                                + result.status
                );
            }
        }
    }

    public static void writeGenes(File file, List<GeneResult> results) throws IOException {
        try (PrintWriter writer = writer(file)) {
            writer.println(
                    "gene_id\tgene_name\tbiotype\tchrom\tstrand\tn_terminal_branches\t"
                            + "control_ale_wpi\tcase_ale_wpi\tdelta_ale_wpi\t"
                            + "control_ale_wpi_sd\tcase_ale_wpi_sd\tale_switch_index\tale_min_pvalue\tale_min_fdr\t"
                            + "control_n_total\tcase_n_total\tcontrol_n_any_evidence\tcase_n_any_evidence\t"
                            + "control_n_informative\tcase_n_informative\tale_change\t"
                            + "observed_ale_direction\tale_direction\tstatus"
            );
            for (GeneResult result : results) {
                writer.println(
                        result.group.geneId + "\t"
                                + result.group.geneName + "\t"
                                + result.group.biotype + "\t"
                                + result.group.chrom + "\t"
                                + result.group.strand.symbol() + "\t"
                                + result.group.branches.size() + "\t"
                                + fmt(result.controlWpi, "%.8f") + "\t"
                                + fmt(result.caseWpi, "%.8f") + "\t"
                                + fmt(result.deltaWpi, "%.8f") + "\t"
                                + fmt(result.controlWpiSd, "%.8f") + "\t"
                                + fmt(result.caseWpiSd, "%.8f") + "\t"
                                + fmt(result.switchIndex, "%.8f") + "\t"
                                + fmt(result.minPvalue, "%.10g") + "\t"
                                + fmt(result.minFdr, "%.10g") + "\t"
                                + result.controlTotalReplicates + "\t"
                                + result.caseTotalReplicates + "\t"
                                + result.controlAnyEvidence + "\t"
                                + result.caseAnyEvidence + "\t"
                                + result.controlInformativeReplicates + "\t"
                                + result.caseInformativeReplicates + "\t"
                                + (result.aleChange ? "yes" : "no") + "\t"
                                + result.observedDirection + "\t"
                                + result.aleDirection + "\t"
                                + result.status
                );
            }
        }
    }

    private static String positions1(List<Integer> positions) {
        StringBuilder out = new StringBuilder();
        for (Integer position : positions) {
            if (out.length() > 0) out.append(',');
            out.append(position + 1);
        }
        return out.toString();
    }

    private static int median(List<Integer> values) {
        List<Integer> copy = new ArrayList<Integer>(values);
        Collections.sort(copy);
        return copy.get(copy.size() / 2);
    }

    private static Map<String, List<TranscriptModel>> groupByGene(
            List<TranscriptModel> transcripts
    ) {
        Map<String, List<TranscriptModel>> grouped =
                new LinkedHashMap<String, List<TranscriptModel>>();
        for (TranscriptModel transcript : transcripts) {
            grouped.computeIfAbsent(
                    transcript.geneId,
                    key -> new ArrayList<TranscriptModel>()
            ).add(transcript);
        }
        return grouped;
    }

    private static String preferredBiotype(
            List<TerminalFamilyReconciler.Family> families
    ) {
        String biotype = ".";
        for (TerminalFamilyReconciler.Family family : families) {
            if ("protein_coding".equals(family.biotype)) return "protein_coding";
            if (".".equals(biotype) && !".".equals(family.biotype)) {
                biotype = family.biotype;
            }
        }
        return biotype;
    }

    private static void addInPlace(double[] target, double[] values) {
        for (int index = 0; index < target.length; index++) target[index] += values[index];
    }

    private static void averageInPlace(double[] values, int count) {
        for (int index = 0; index < values.length; index++) values[index] /= count;
    }

    private static double totalVariation(double[] left, double[] right) {
        double sum = 0.0;
        for (int index = 0; index < left.length; index++) {
            sum += Math.abs(right[index] - left[index]);
        }
        return 0.5 * sum;
    }

    private static double mean(List<Double> values) {
        if (values.isEmpty()) return Double.NaN;
        double sum = 0.0;
        for (Double value : values) sum += value;
        return sum / values.size();
    }

    private static double sampleSd(List<Double> values) {
        if (values.size() < 2) return Double.NaN;
        double mean = mean(values);
        double ss = 0.0;
        for (Double value : values) {
            double delta = value - mean;
            ss += delta * delta;
        }
        return Math.sqrt(ss / (values.size() - 1.0));
    }

    private static double finiteDifference(double right, double left) {
        return Double.isFinite(right) && Double.isFinite(left)
                ? right - left
                : Double.NaN;
    }

    private static int compareFinite(double left, double right) {
        boolean lf = Double.isFinite(left);
        boolean rf = Double.isFinite(right);
        if (lf && rf) return Double.compare(left, right);
        if (lf) return -1;
        if (rf) return 1;
        return 0;
    }

    private static String fmt(double value, String format) {
        return Double.isFinite(value)
                ? String.format(Locale.ROOT, format, value)
                : "NA";
    }

    private static PrintWriter writer(File file) throws IOException {
        return new PrintWriter(new BufferedWriter(new FileWriter(file)));
    }
}
