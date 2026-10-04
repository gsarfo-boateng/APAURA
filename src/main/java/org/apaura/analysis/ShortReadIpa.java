package org.apaura.analysis;

import org.apaura.model.CoverageSegment;
import org.apaura.model.Interval;
import org.apaura.model.Sample;
import org.apaura.model.Strand;
import org.apaura.model.TranscriptModel;
import org.apaura.stats.BetaBinomial;
import org.apaura.stats.FisherExact;
import org.apaura.stats.MultipleTesting;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/**
 * Annotated intronic polyadenylation analysis for conventional short-read RNA-seq.
 *
 * <p>An IPA candidate is an annotated transcript end that falls inside an intron of
 * another transcript from the same gene. Quantification uses equal-length exonic
 * windows around that host intron: a shared upstream exon window provides QC,
 * an intronic terminal-exon window is specific to the annotated IPA isoform, and a
 * downstream exon window represents transcripts that continue beyond the intron.
 * IPA usage is IPA-specific / (IPA-specific + continuation).</p>
 */
public final class ShortReadIpa {

    private static final int DEFAULT_MIN_INFORMATIVE_REPLICATES = 3;
    private ShortReadIpa() {
        // Utility class.
    }

    public static final class Group {
        public final String id;
        public final String geneId;
        public final String geneName;
        public final String biotype;
        public final String chrom;
        public final Strand strand;
        public final int intronStart0;
        public final int intronEnd0;
        public final int representativePas0;
        public final List<Integer> annotatedPas0;
        public final CoverageSegment commonWindow;
        public final CoverageSegment ipaSpecificWindow;
        public final CoverageSegment continuationWindow;

        Group(
                String id,
                String geneId,
                String geneName,
                String biotype,
                String chrom,
                Strand strand,
                int intronStart0,
                int intronEnd0,
                int representativePas0,
                List<Integer> annotatedPas0,
                CoverageSegment commonWindow,
                CoverageSegment ipaSpecificWindow,
                CoverageSegment continuationWindow
        ) {
            this.id = id;
            this.geneId = geneId;
            this.geneName = geneName;
            this.biotype = biotype;
            this.chrom = chrom;
            this.strand = strand;
            this.intronStart0 = intronStart0;
            this.intronEnd0 = intronEnd0;
            this.representativePas0 = representativePas0;
            this.annotatedPas0 = Collections.unmodifiableList(
                    new ArrayList<Integer>(annotatedPas0)
            );
            this.commonWindow = commonWindow;
            this.ipaSpecificWindow = ipaSpecificWindow;
            this.continuationWindow = continuationWindow;
        }
    }

    public static final class Quantification {
        public final Group group;
        public final int commonReads;
        public final int ipaSpecificReads;
        public final int continuationReads;
        public final int totalEvidenceReads;
        public final double ipaUsage;

        Quantification(
                Group group,
                int commonReads,
                int ipaSpecificReads,
                int continuationReads
        ) {
            this.group = group;
            this.commonReads = commonReads;
            this.ipaSpecificReads = ipaSpecificReads;
            this.continuationReads = continuationReads;
            this.totalEvidenceReads = ipaSpecificReads + continuationReads;
            this.ipaUsage = totalEvidenceReads > 0
                    ? (double) ipaSpecificReads / totalEvidenceReads
                    : 0.0;
        }
    }

    public static final class SiteResult {
        public Group group;
        public double controlUsage = Double.NaN;
        public double caseUsage = Double.NaN;
        public double delta = Double.NaN;
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
        public String direction = "not_applicable";
        public String method = "not_tested";
        public String status = "insufficient_evidence";
    }

    public static final class GeneResult {
        public String geneId;
        public String geneName;
        public String biotype;
        public String chrom;
        public Strand strand;
        public int nIpaEvents;
        public int nSignificantEvents;
        public String ipaDirection = "stable";
        public double maxAbsDelta;
        public double minFdr = 1.0;
    }

    public static final class Results {
        public final List<SiteResult> sites = new ArrayList<SiteResult>();
        public final List<GeneResult> genes = new ArrayList<GeneResult>();
    }

    private static final class Support {
        TranscriptModel candidate;
        int pas0;
        int intronStart0;
        int intronEnd0;
        CoverageSegment commonWindow;
        CoverageSegment ipaSpecificWindow;
        CoverageSegment continuationWindow;

        int windowLength() {
            return commonWindow.length();
        }

        int intronLength() {
            return intronEnd0 - intronStart0;
        }
    }

    private static final class EventBuilder {
        final Support support;
        final Set<Integer> pasPositions = new LinkedHashSet<Integer>();

        EventBuilder(Support support) {
            this.support = support;
        }
    }

    /** Build non-redundant, annotation-supported IPA events. */
    public static List<Group> build(
            List<TranscriptModel> transcripts,
            int pasMergeDistance,
            int windowSize,
            int minWindowLength
    ) {
        if (pasMergeDistance < 0) {
            throw new IllegalArgumentException(
                    "--pas-merge-distance must be >= 0"
            );
        }
        if (windowSize < minWindowLength) {
            throw new IllegalArgumentException(
                    "--ipa-window must be >= --min-segment-length"
            );
        }

        Map<String, List<TranscriptModel>> byGene = groupByGene(transcripts);
        List<Group> groups = new ArrayList<Group>();

        for (List<TranscriptModel> geneTranscripts : byGene.values()) {
            List<Group> geneGroups = buildGene(
                    geneTranscripts,
                    pasMergeDistance,
                    windowSize,
                    minWindowLength
            );
            groups.addAll(geneGroups);
        }

        groups.sort(new Comparator<Group>() {
            @Override
            public int compare(Group left, Group right) {
                int comparison = left.chrom.compareTo(right.chrom);
                if (comparison != 0) {
                    return comparison;
                }
                comparison = left.geneId.compareTo(right.geneId);
                if (comparison != 0) {
                    return comparison;
                }
                return Integer.compare(left.representativePas0, right.representativePas0);
            }
        });
        return Collections.unmodifiableList(groups);
    }

    private static List<Group> buildGene(
            List<TranscriptModel> transcripts,
            int pasMergeDistance,
            int windowSize,
            int minWindowLength
    ) {
        if (transcripts.isEmpty()) {
            return Collections.emptyList();
        }

        Map<String, Support> bestSupportByPas = new LinkedHashMap<String, Support>();
        for (TranscriptModel candidate : transcripts) {
            int pas0 = candidate.pasPos0();
            Support best = null;

            for (TranscriptModel host : transcripts) {
                if (candidate == host
                        || !candidate.chrom.equals(host.chrom)
                        || candidate.strand != host.strand) {
                    continue;
                }

                Support support = support(
                        candidate,
                        host,
                        pas0,
                        windowSize,
                        minWindowLength
                );
                if (support != null && better(support, best)) {
                    best = support;
                }
            }

            if (best != null) {
                String key = candidate.geneId
                        + "\u0001" + candidate.chrom
                        + "\u0001" + candidate.strand.symbol()
                        + "\u0001" + pas0;
                Support previous = bestSupportByPas.get(key);
                if (better(best, previous)) {
                    bestSupportByPas.put(key, best);
                }
            }
        }

        /*
         * A host intron can contain several distinct annotated cleavage loci.
         * Keep the intron context, but split its supported transcript ends into
         * local physical PAS clusters using the same --pas-merge-distance used
         * by terminal APA. This prevents distant IPA sites in one intron from
         * being collapsed into a single event.
         */
        Map<String, List<Support>> supportsByHostIntron =
                new LinkedHashMap<String, List<Support>>();
        for (Support support : bestSupportByPas.values()) {
            String key = support.candidate.geneId
                    + "\u0001" + support.candidate.chrom
                    + "\u0001" + support.candidate.strand.symbol()
                    + "\u0001" + support.intronStart0
                    + "\u0001" + support.intronEnd0;
            supportsByHostIntron.computeIfAbsent(
                    key,
                    ignored -> new ArrayList<Support>()
            ).add(support);
        }

        List<EventBuilder> events = new ArrayList<EventBuilder>();
        for (List<Support> intronSupports : supportsByHostIntron.values()) {
            for (List<Support> cluster : clusterSupports(
                    intronSupports,
                    pasMergeDistance
            )) {
                int representativePas0 = representativePas0(cluster);
                Support representativeSupport = supportAtPas(
                        cluster,
                        representativePas0
                );
                EventBuilder event = new EventBuilder(representativeSupport);
                for (Support support : cluster) {
                    event.pasPositions.add(support.pas0);
                }
                events.add(event);
            }
        }

        events.sort(new Comparator<EventBuilder>() {
            @Override
            public int compare(EventBuilder left, EventBuilder right) {
                int comparison = Integer.compare(
                        left.support.intronStart0,
                        right.support.intronStart0
                );
                if (comparison != 0) {
                    return comparison;
                }
                comparison = Integer.compare(
                        left.support.intronEnd0,
                        right.support.intronEnd0
                );
                if (comparison != 0) {
                    return comparison;
                }
                return Integer.compare(
                        representativePas0(left),
                        representativePas0(right)
                );
            }
        });

        TranscriptModel first = transcripts.get(0);
        String biotype = preferredBiotype(transcripts);
        List<Group> groups = new ArrayList<Group>();
        for (int index = 0; index < events.size(); index++) {
            EventBuilder event = events.get(index);
            List<Integer> positions = new ArrayList<Integer>(event.pasPositions);
            Collections.sort(positions);
            int representative = representativePas0(event);
            String id = first.geneId
                    + "|IPA"
                    + String.format(Locale.ROOT, "%03d", index + 1);

            groups.add(new Group(
                    id,
                    first.geneId,
                    first.geneName,
                    biotype,
                    first.chrom,
                    first.strand,
                    event.support.intronStart0,
                    event.support.intronEnd0,
                    representative,
                    positions,
                    event.support.commonWindow,
                    event.support.ipaSpecificWindow,
                    event.support.continuationWindow
            ));
        }
        return groups;
    }

    private static List<List<Support>> clusterSupports(
            List<Support> supports,
            int mergeDistance
    ) {
        List<Support> ordered = new ArrayList<Support>(supports);
        ordered.sort(new Comparator<Support>() {
            @Override
            public int compare(Support left, Support right) {
                return Integer.compare(left.pas0, right.pas0);
            }
        });

        List<List<Support>> clusters = new ArrayList<List<Support>>();
        if (ordered.isEmpty()) {
            return clusters;
        }

        List<Support> cluster = new ArrayList<Support>();
        cluster.add(ordered.get(0));
        for (int index = 1; index < ordered.size(); index++) {
            Support current = ordered.get(index);
            Support previous = ordered.get(index - 1);
            if (current.pas0 - previous.pas0 <= mergeDistance) {
                cluster.add(current);
            } else {
                clusters.add(cluster);
                cluster = new ArrayList<Support>();
                cluster.add(current);
            }
        }
        clusters.add(cluster);
        return clusters;
    }

    private static int representativePas0(List<Support> cluster) {
        List<Integer> positions = new ArrayList<Integer>();
        for (Support support : cluster) {
            positions.add(support.pas0);
        }
        Collections.sort(positions);
        return positions.get(positions.size() / 2);
    }

    private static int representativePas0(EventBuilder event) {
        List<Integer> positions = new ArrayList<Integer>(event.pasPositions);
        Collections.sort(positions);
        return positions.get(positions.size() / 2);
    }

    private static Support supportAtPas(
            List<Support> supports,
            int pas0
    ) {
        Support best = null;
        for (Support support : supports) {
            if (support.pas0 == pas0 && better(support, best)) {
                best = support;
            }
        }
        if (best != null) {
            return best;
        }

        for (Support support : supports) {
            if (better(support, best)) {
                best = support;
            }
        }
        return best;
    }

    private static Support support(
            TranscriptModel candidate,
            TranscriptModel host,
            int pas0,
            int windowSize,
            int minWindowLength
    ) {
        List<Interval> exons = sortedExons(host);
        for (int index = 1; index < exons.size(); index++) {
            Interval left = exons.get(index - 1);
            Interval right = exons.get(index);
            if (left.end0 >= right.start0) {
                continue;
            }
            int intronStart0 = left.end0;
            int intronEnd0 = right.start0;
            if (pas0 < intronStart0 || pas0 >= intronEnd0) {
                continue;
            }

            Interval upstream = host.strand == Strand.MINUS ? right : left;
            Interval continuation = host.strand == Strand.MINUS ? left : right;
            Interval shared = bestSharedExonicInterval(candidate, upstream, host.strand);
            Interval terminal = candidate.terminalExon();
            Interval ipaSpecific = intersect(terminal, new Interval(intronStart0, intronEnd0));
            if (shared == null || ipaSpecific == null) {
                continue;
            }

            int length = Math.min(
                    windowSize,
                    Math.min(
                            Math.min(
                                    shared.end0 - shared.start0,
                                    ipaSpecific.end0 - ipaSpecific.start0
                            ),
                            continuation.end0 - continuation.start0
                    )
            );
            if (length < minWindowLength) {
                continue;
            }

            CoverageSegment commonWindow;
            CoverageSegment ipaSpecificWindow;
            CoverageSegment continuationWindow;
            if (host.strand == Strand.MINUS) {
                commonWindow = new CoverageSegment(
                        0,
                        shared.start0,
                        shared.start0 + length,
                        "COMMON_UPSTREAM"
                );
                ipaSpecificWindow = new CoverageSegment(
                        1,
                        ipaSpecific.start0,
                        ipaSpecific.start0 + length,
                        "IPA_SPECIFIC"
                );
                continuationWindow = new CoverageSegment(
                        2,
                        continuation.end0 - length,
                        continuation.end0,
                        "CONTINUATION"
                );
            } else {
                commonWindow = new CoverageSegment(
                        0,
                        shared.end0 - length,
                        shared.end0,
                        "COMMON_UPSTREAM"
                );
                ipaSpecificWindow = new CoverageSegment(
                        1,
                        ipaSpecific.end0 - length,
                        ipaSpecific.end0,
                        "IPA_SPECIFIC"
                );
                continuationWindow = new CoverageSegment(
                        2,
                        continuation.start0,
                        continuation.start0 + length,
                        "CONTINUATION"
                );
            }

            Support support = new Support();
            support.candidate = candidate;
            support.pas0 = pas0;
            support.intronStart0 = intronStart0;
            support.intronEnd0 = intronEnd0;
            support.commonWindow = commonWindow;
            support.ipaSpecificWindow = ipaSpecificWindow;
            support.continuationWindow = continuationWindow;
            return support;
        }
        return null;
    }

    private static Interval intersect(Interval left, Interval right) {
        if (left == null || right == null) {
            return null;
        }
        int start0 = Math.max(left.start0, right.start0);
        int end0 = Math.min(left.end0, right.end0);
        return end0 > start0 ? new Interval(start0, end0) : null;
    }

    private static Interval bestSharedExonicInterval(
            TranscriptModel transcript,
            Interval hostUpstreamExon,
            Strand strand
    ) {
        Interval best = null;
        for (Interval exon : transcript.exons) {
            int start0 = Math.max(exon.start0, hostUpstreamExon.start0);
            int end0 = Math.min(exon.end0, hostUpstreamExon.end0);
            if (end0 <= start0) {
                continue;
            }
            Interval overlap = new Interval(start0, end0);
            if (best == null || betterShared(overlap, best, strand)) {
                best = overlap;
            }
        }
        return best;
    }

    private static boolean betterShared(Interval candidate, Interval current, Strand strand) {
        if ((candidate.end0 - candidate.start0) != (current.end0 - current.start0)) {
            return (candidate.end0 - candidate.start0) > (current.end0 - current.start0);
        }
        return strand == Strand.MINUS
                ? candidate.start0 < current.start0
                : candidate.end0 > current.end0;
    }

    private static boolean better(Support candidate, Support current) {
        if (candidate == null) {
            return false;
        }
        if (current == null) {
            return true;
        }
        if (candidate.windowLength() != current.windowLength()) {
            return candidate.windowLength() > current.windowLength();
        }
        return candidate.intronLength() < current.intronLength();
    }

    private static List<Interval> sortedExons(TranscriptModel transcript) {
        List<Interval> exons = new ArrayList<Interval>(transcript.exons);
        exons.sort(Comparator.comparingInt(interval -> interval.start0));
        return exons;
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

    private static String preferredBiotype(List<TranscriptModel> transcripts) {
        String biotype = ".";
        for (TranscriptModel transcript : transcripts) {
            if ("protein_coding".equals(transcript.biotype)) {
                return "protein_coding";
            }
            if (".".equals(biotype) && !".".equals(transcript.biotype)) {
                biotype = transcript.biotype;
            }
        }
        return biotype;
    }

    public static Quantification quantify(
            Group group,
            int commonReads,
            int ipaSpecificReads,
            int continuationReads
    ) {
        return new Quantification(
                group,
                commonReads,
                ipaSpecificReads,
                continuationReads
        );
    }

    public static Results differential(
            List<Sample> samples,
            Map<String, List<Quantification>> bySample,
            String control,
            String caseCondition,
            int minSampleEvidence,
            double minDelta,
            double fdr
    ) {
        return differential(samples, bySample, control, caseCondition, minSampleEvidence,
                DEFAULT_MIN_INFORMATIVE_REPLICATES, minDelta, fdr);
    }

    public static Results differential(
            List<Sample> samples,
            Map<String, List<Quantification>> bySample,
            String control,
            String caseCondition,
            int minSampleEvidence,
            int minInformativeReplicates,
            double minDelta,
            double fdr
    ) {
        Results results = new Results();
        if (samples.isEmpty()) return results;
        List<Quantification> first = bySample.get(samples.get(0).name);
        if (first == null) return results;
        for (int groupIndex = 0; groupIndex < first.size(); groupIndex++) {
            results.sites.add(testEvent(samples, bySample, groupIndex, control, caseCondition,
                    minSampleEvidence, minInformativeReplicates));
        }
        applyFdr(results.sites);
        classify(results.sites, minDelta, fdr);
        summarizeGenes(results);
        return results;
    }

    private static SiteResult testEvent(
            List<Sample> samples,
            Map<String, List<Quantification>> bySample,
            int groupIndex,
            String control,
            String caseCondition,
            int minSampleEvidence,
            int minInformativeReplicates
    ) {
        Group group = bySample.get(samples.get(0).name).get(groupIndex).group;
        int[] ipaCounts = new int[samples.size()];
        int[] totals = new int[samples.size()];
        int[] labels = new int[samples.size()];
        for (int i = 0; i < labels.length; i++) labels[i] = -1;

        double controlSum = 0.0, caseSum = 0.0;
        int controlTotal = 0, caseTotal = 0;
        int controlAny = 0, caseAny = 0;
        int controlInformative = 0, caseInformative = 0;
        int controlEvidence = 0, caseEvidence = 0;

        for (int index = 0; index < samples.size(); index++) {
            Sample sample = samples.get(index);
            Quantification q = bySample.get(sample.name).get(groupIndex);
            boolean isControl = sample.condition.equals(control);
            boolean isCase = sample.condition.equals(caseCondition);
            if (!isControl && !isCase) continue;
            if (isControl) { controlTotal++; if (q.totalEvidenceReads > 0) controlAny++; }
            else { caseTotal++; if (q.totalEvidenceReads > 0) caseAny++; }
            if (q.totalEvidenceReads < minSampleEvidence) continue;
            ipaCounts[index] = q.ipaSpecificReads;
            totals[index] = q.totalEvidenceReads;
            if (isControl) {
                labels[index] = 0; controlSum += q.ipaUsage; controlInformative++; controlEvidence += q.totalEvidenceReads;
            } else {
                labels[index] = 1; caseSum += q.ipaUsage; caseInformative++; caseEvidence += q.totalEvidenceReads;
            }
        }

        SiteResult result = new SiteResult();
        result.group = group;
        result.controlTotalReplicates = controlTotal;
        result.caseTotalReplicates = caseTotal;
        result.controlAnyEvidence = controlAny;
        result.caseAnyEvidence = caseAny;
        result.controlInformativeReplicates = controlInformative;
        result.caseInformativeReplicates = caseInformative;
        result.controlEvidence = controlEvidence;
        result.caseEvidence = caseEvidence;
        if (controlInformative == 0 || caseInformative == 0) return result;

        result.controlUsage = controlSum / controlInformative;
        result.caseUsage = caseSum / caseInformative;
        result.delta = result.caseUsage - result.controlUsage;

        boolean oneVsOne = controlTotal == 1 && caseTotal == 1;
        if (oneVsOne) {
            int[] table = pooledTwoByTwo(ipaCounts, totals, labels);
            result.p = FisherExact.twoSided(table[0], table[1], table[2], table[3]);
            result.method = "IPA-coverage-Fisher-exact";
            result.status = "ok";
            return result;
        }
        if (controlInformative < minInformativeReplicates || caseInformative < minInformativeReplicates) {
            result.status = "insufficient_replicate_support";
            return result;
        }
        BetaBinomial.TestResult test = BetaBinomial.test(ipaCounts, totals, labels);
        result.p = test.pvalue;
        result.rho = test.rho;
        result.method = "IPA-coverage-beta-binomial-LRT";
        result.status = "ok";
        return result;
    }

    private static int[] pooledTwoByTwo(int[] successes, int[] totals, int[] labels) {
        int controlSuccess = 0;
        int controlOther = 0;
        int caseSuccess = 0;
        int caseOther = 0;
        for (int index = 0; index < labels.length; index++) {
            if (labels[index] == 0) {
                controlSuccess += successes[index];
                controlOther += totals[index] - successes[index];
            } else if (labels[index] == 1) {
                caseSuccess += successes[index];
                caseOther += totals[index] - successes[index];
            }
        }
        return new int[]{controlSuccess, controlOther, caseSuccess, caseOther};
    }

    private static void applyFdr(List<SiteResult> results) {
        List<Integer> idx = new ArrayList<Integer>();
        List<Double> ps = new ArrayList<Double>();
        for (int i = 0; i < results.size(); i++) {
            if (Double.isFinite(results.get(i).p)) { idx.add(i); ps.add(results.get(i).p); }
        }
        if (ps.isEmpty()) return;
        double[] pValues = new double[ps.size()];
        for (int i = 0; i < ps.size(); i++) pValues[i] = ps.get(i);
        double[] qValues = MultipleTesting.bh(pValues);
        for (int i = 0; i < qValues.length; i++) results.get(idx.get(i)).q = qValues[i];
    }

    private static void classify(List<SiteResult> results, double minDelta, double fdr) {
        for (SiteResult result : results) {
            if (!"ok".equals(result.status)) {
                result.direction = "not_applicable";
            } else if (result.q <= fdr && Math.abs(result.delta) >= minDelta) {
                result.direction = result.delta > 0.0 ? "gain" : "loss";
            } else {
                result.direction = "stable";
            }
        }
    }

    private static void summarizeGenes(Results results) {
        Map<String, GeneResult> byGene = new LinkedHashMap<String, GeneResult>();
        Map<String, Boolean> hasGain = new LinkedHashMap<String, Boolean>();
        Map<String, Boolean> hasLoss = new LinkedHashMap<String, Boolean>();

        for (SiteResult site : results.sites) {
            GeneResult gene = byGene.get(site.group.geneId);
            if (gene == null) {
                gene = new GeneResult();
                gene.geneId = site.group.geneId;
                gene.geneName = site.group.geneName;
                gene.biotype = site.group.biotype;
                gene.chrom = site.group.chrom;
                gene.strand = site.group.strand;
                byGene.put(gene.geneId, gene);
            }
            gene.nIpaEvents++;
            gene.minFdr = Math.min(gene.minFdr, site.q);
            gene.maxAbsDelta = Math.max(gene.maxAbsDelta, Math.abs(site.delta));
            if ("gain".equals(site.direction)) {
                gene.nSignificantEvents++;
                hasGain.put(gene.geneId, true);
            } else if ("loss".equals(site.direction)) {
                gene.nSignificantEvents++;
                hasLoss.put(gene.geneId, true);
            }
        }

        for (GeneResult gene : byGene.values()) {
            boolean gain = Boolean.TRUE.equals(hasGain.get(gene.geneId));
            boolean loss = Boolean.TRUE.equals(hasLoss.get(gene.geneId));
            if (gain && loss) {
                gene.ipaDirection = "complex";
            } else if (gain) {
                gene.ipaDirection = "gain";
            } else if (loss) {
                gene.ipaDirection = "loss";
            }
            results.genes.add(gene);
        }
    }

    public static void writeAnnotation(File file, List<Group> groups) throws IOException {
        try (PrintWriter writer = writer(file)) {
            writer.println(
                    "ipa_id\tgene_id\tgene_name\tbiotype\tchrom\tstrand\tintron_start0\t"
                            + "intron_end0\trepresentative_pas_1based\tn_annotated_ipa_pas\t"
                            + "annotated_ipa_positions_1based\tcommon_start0\tcommon_end0\t"
                            + "ipa_specific_start0\tipa_specific_end0\tcontinuation_start0\t"
                            + "continuation_end0\twindow_length"
            );
            for (Group group : groups) {
                writer.println(annotationLine(group));
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
                    "sample\tcondition\tipa_id\tgene_id\tgene_name\tchrom\tstrand\t"
                            + "representative_pas_1based\tcommon_reads\tipa_specific_reads\t"
                            + "continuation_reads\tipa_usage\ttotal_evidence_reads"
            );
            for (Sample sample : samples) {
                for (Quantification quantification : bySample.get(sample.name)) {
                    writer.printf(
                            Locale.ROOT,
                            "%s\t%s\t%s\t%s\t%s\t%s\t%c\t%d\t%d\t%d\t%d\t%.8f\t%d%n",
                            sample.name,
                            sample.condition,
                            quantification.group.id,
                            quantification.group.geneId,
                            quantification.group.geneName,
                            quantification.group.chrom,
                            quantification.group.strand.symbol(),
                            quantification.group.representativePas0 + 1,
                            quantification.commonReads,
                            quantification.ipaSpecificReads,
                            quantification.continuationReads,
                            quantification.ipaUsage,
                            quantification.totalEvidenceReads
                    );
                }
            }
        }
    }

    public static void writeSites(File file, List<SiteResult> results) throws IOException {
        try (PrintWriter writer = writer(file)) {
            writer.println(
                    "ipa_id\tgene_id\tgene_name\tbiotype\tchrom\tstrand\tintron_start0\t"
                            + "intron_end0\trepresentative_pas_1based\tn_annotated_ipa_pas\t"
                            + "annotated_ipa_positions_1based\tcontrol_ipa_usage\tcase_ipa_usage\t"
                            + "delta_ipa_usage\tcontrol_n_total\tcase_n_total\t"
                            + "control_n_any_evidence\tcase_n_any_evidence\t"
                            + "control_n_informative\tcase_n_informative\t"
                            + "control_evidence\tcase_evidence\tpvalue\tfdr\t"
                            + "dispersion_rho\tipa_direction\tmethod\tstatus"
            );

            for (SiteResult result : results) {
                Group group = result.group;

                writer.println(
                        group.id + "\t"
                                + group.geneId + "\t"
                                + group.geneName + "\t"
                                + group.biotype + "\t"
                                + group.chrom + "\t"
                                + group.strand.symbol() + "\t"
                                + group.intronStart0 + "\t"
                                + group.intronEnd0 + "\t"
                                + (group.representativePas0 + 1) + "\t"
                                + group.annotatedPas0.size() + "\t"
                                + positions(group.annotatedPas0) + "\t"
                                + fmt(result.controlUsage, "%.8f") + "\t"
                                + fmt(result.caseUsage, "%.8f") + "\t"
                                + fmt(result.delta, "%.8f") + "\t"
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
                                + result.direction + "\t"
                                + result.method + "\t"
                                + result.status
                );
            }
        }
    }

    private static String fmt(double value, String format) {
        return Double.isFinite(value)
                ? String.format(Locale.ROOT, format, value)
                : "NA";
    }

    public static void writeGenes(File file, List<GeneResult> results) throws IOException {
        try (PrintWriter writer = writer(file)) {
            writer.println(
                    "gene_id\tgene_name\tbiotype\tchrom\tstrand\tn_ipa_events\t"
                            + "n_significant_ipa_events\tipa_direction\tmax_abs_delta_ipa\tmin_fdr"
            );
            for (GeneResult result : results) {
                writer.printf(
                        Locale.ROOT,
                        "%s\t%s\t%s\t%s\t%c\t%d\t%d\t%s\t%.8f\t%.10g%n",
                        result.geneId,
                        result.geneName,
                        result.biotype,
                        result.chrom,
                        result.strand.symbol(),
                        result.nIpaEvents,
                        result.nSignificantEvents,
                        result.ipaDirection,
                        result.maxAbsDelta,
                        result.minFdr
                );
            }
        }
    }

    private static String annotationLine(Group group) {
        return group.id + "\t"
                + group.geneId + "\t"
                + group.geneName + "\t"
                + group.biotype + "\t"
                + group.chrom + "\t"
                + group.strand.symbol() + "\t"
                + group.intronStart0 + "\t"
                + group.intronEnd0 + "\t"
                + (group.representativePas0 + 1) + "\t"
                + group.annotatedPas0.size() + "\t"
                + positions(group.annotatedPas0) + "\t"
                + group.commonWindow.start0 + "\t"
                + group.commonWindow.end0 + "\t"
                + group.ipaSpecificWindow.start0 + "\t"
                + group.ipaSpecificWindow.end0 + "\t"
                + group.continuationWindow.start0 + "\t"
                + group.continuationWindow.end0 + "\t"
                + group.commonWindow.length();
    }

    private static String positions(List<Integer> positions0) {
        StringBuilder values = new StringBuilder();
        for (Integer position0 : positions0) {
            if (values.length() > 0) {
                values.append(',');
            }
            values.append(position0 + 1);
        }
        return values.toString();
    }

    private static PrintWriter writer(File file) throws IOException {
        return new PrintWriter(new BufferedWriter(new FileWriter(file)));
    }
}