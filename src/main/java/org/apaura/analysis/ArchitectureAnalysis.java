package org.apaura.analysis;

import org.apaura.model.Sample;
import org.apaura.model.SiteCluster;
import org.apaura.model.Strand;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Computes gene-level APA architecture metrics from annotated PAS counts.
 * The analysis separates termination compartment usage, alternative last-exon
 * branch choice, and PAS choice within terminal-exon branches.
 */
public final class ArchitectureAnalysis {
    private static final int TERMINAL = 0;
    private static final int IPA = 1;
    private static final int INTERNAL = 2;
    private static final int DOWNSTREAM = 3;
    private static final int COMPARTMENT_COUNT = 4;

    private static final String[] COMPARTMENTS = {
            "terminal",
            "IPA",
            "internal_exon",
            "downstream"
    };

    private ArchitectureAnalysis() {
    }

    public static final class Metrics {
        public String geneId;
        public String architectureClass;
        public String components;
        public String dominantControl;
        public String dominantCase;
        public String transition;
        public String mechanismEvents;
        public int nTerminalBranches;
        public int nTutrBranches;
        public double controlTerminalUsage;
        public double caseTerminalUsage;
        public double deltaTerminalUsage;
        public double controlIpaUsage;
        public double caseIpaUsage;
        public double deltaIpaUsage;
        public double controlInternalUsage;
        public double caseInternalUsage;
        public double deltaInternalUsage;
        public double controlDownstreamUsage;
        public double caseDownstreamUsage;
        public double deltaDownstreamUsage;
        public double terminationShiftIndex;
        public double terminationJsd;
        public double aleSwitchIndex;
        public double aleJsd;
        public double controlAleWpi;
        public double caseAleWpi;
        public double deltaAleWpi;
        public double controlTutrWpi;
        public double caseTutrWpi;
        public double deltaTutrWpi;
        public double tutrShiftMagnitude;
    }

    private static final class GeneAccumulator {
        final Map<Integer, List<SiteCluster>> branches;
        final List<Integer> anchors;
        final Strand strand;

        final double[] controlCompartments = new double[COMPARTMENT_COUNT];
        final double[] caseCompartments = new double[COMPARTMENT_COUNT];
        int controlCompartmentSamples;
        int caseCompartmentSamples;

        final double[] controlBranches;
        final double[] caseBranches;
        int controlBranchSamples;
        int caseBranchSamples;

        final double[] controlTutrWpi;
        final double[] caseTutrWpi;
        final int[] controlTutrSamples;
        final int[] caseTutrSamples;
        final List<double[]> controlTutrPau;
        final List<double[]> caseTutrPau;

        GeneAccumulator(
                Map<Integer, List<SiteCluster>> branches,
                List<Integer> anchors,
                Strand strand
        ) {
            this.branches = branches;
            this.anchors = anchors;
            this.strand = strand;
            this.controlBranches = new double[anchors.size()];
            this.caseBranches = new double[anchors.size()];
            this.controlTutrWpi = new double[anchors.size()];
            this.caseTutrWpi = new double[anchors.size()];
            this.controlTutrSamples = new int[anchors.size()];
            this.caseTutrSamples = new int[anchors.size()];
            this.controlTutrPau = new ArrayList<double[]>();
            this.caseTutrPau = new ArrayList<double[]>();

            for (Integer anchor : anchors) {
                int siteCount = branches.get(anchor).size();
                controlTutrPau.add(new double[siteCount]);
                caseTutrPau.add(new double[siteCount]);
            }
        }
    }

    public static void annotateSiteRoles(List<SiteCluster> sites) {
        Map<String, List<SiteCluster>> sitesByGene = groupByGene(sites);

        for (List<SiteCluster> geneSites : sitesByGene.values()) {
            Strand strand = geneSites.get(0).strand;
            Map<Integer, List<SiteCluster>> branches = terminalBranches(geneSites);
            List<Integer> anchors = new ArrayList<Integer>(branches.keySet());
            Collections.sort(anchors);
            if (strand == Strand.MINUS) {
                Collections.reverse(anchors);
            }

            Map<Integer, Integer> branchRanks = new HashMap<Integer, Integer>();
            for (int i = 0; i < anchors.size(); i++) {
                branchRanks.put(anchors.get(i), i + 1);
            }

            boolean hasAle = anchors.size() > 1;
            for (SiteCluster site : geneSites) {
                assignSiteRole(site, branches, branchRanks, hasAle);
            }
        }
    }

    private static void assignSiteRole(
            SiteCluster site,
            Map<Integer, List<SiteCluster>> branches,
            Map<Integer, Integer> branchRanks,
            boolean hasAle
    ) {
        if (site.terminalExon &&
                site.terminalStructureAnchor0 >= 0) {
            Integer rank = branchRanks.get(
                    site.terminalStructureAnchor0
            );
            site.terminalBranchRank =
                    rank == null ? -1 : rank;

            int branchSize = branches.get(
                    site.terminalStructureAnchor0
            ).size();

            if (hasAle && branchSize > 1) {
                site.architectureRole = "ALE+TUTR";
            } else if (hasAle) {
                site.architectureRole = "ALE";
            } else if (branchSize > 1) {
                site.architectureRole = "TUTR";
            } else {
                site.architectureRole = "terminal";
            }
            return;
        }

        if ("INTRON".equals(site.region)) {
            site.architectureRole = "IPA";
        } else if ("EXON_CODING".equals(site.region)
                || "EXON_NONCODING".equals(site.region)) {
            site.architectureRole = "internal_exon";
        } else if ("DOWNSTREAM_3P".equals(site.region)) {
            site.architectureRole = "downstream";
        } else {
            site.architectureRole = "other";
        }
    }

    public static Map<String, Metrics> analyze(
            List<Sample> samples,
            List<SiteCluster> sites,
            String control,
            String caseCondition
    ) {
        annotateSiteRoles(sites);
        Map<String, Metrics> results = new LinkedHashMap<String, Metrics>();

        for (Map.Entry<String, List<SiteCluster>> entry : groupByGene(sites).entrySet()) {
            if (entry.getValue().isEmpty()) {
                continue;
            }
            results.put(
                    entry.getKey(),
                    analyzeGene(
                            entry.getKey(),
                            entry.getValue(),
                            samples,
                            control,
                            caseCondition
                    )
            );
        }
        return results;
    }

    private static Metrics analyzeGene(
            String geneId,
            List<SiteCluster> geneSites,
            List<Sample> samples,
            String control,
            String caseCondition
    ) {
        Metrics metrics = new Metrics();
        metrics.geneId = geneId;

        Map<Integer, List<SiteCluster>> branches = terminalBranches(geneSites);
        metrics.nTerminalBranches = branches.size();
        metrics.nTutrBranches = countTutrBranches(branches);

        boolean hasAle = branches.size() >= 2;
        boolean hasTutr = metrics.nTutrBranches > 0;
        metrics.components = architectureComponents(hasTutr, hasAle);
        metrics.architectureClass = hasUnresolvedTerminalPas(geneSites)
                ? ApaArchitecture.UNRESOLVED
                : architectureClass(
                        metrics.nTerminalBranches,
                        metrics.nTutrBranches,
                        countTerminalPas(branches)
                );

        Strand strand = geneSites.get(0).strand;
        List<Integer> anchors = orderedAnchors(branches, strand);
        GeneAccumulator accumulator = new GeneAccumulator(branches, anchors, strand);

        for (Sample sample : samples) {
            accumulateSample(
                    accumulator,
                    geneSites,
                    sample,
                    control,
                    caseCondition
            );
        }

        finalizeMetrics(metrics, accumulator);
        return metrics;
    }

    private static int countTutrBranches(Map<Integer, List<SiteCluster>> branches) {
        int count = 0;
        for (List<SiteCluster> branch : branches.values()) {
            if (branch.size() >= 2) {
                count++;
            }
        }
        return count;
    }

    private static int countTerminalPas(
            Map<Integer, List<SiteCluster>> branches
    ) {
        int count = 0;
        for (List<SiteCluster> sites : branches.values()) {
            count += sites.size();
        }
        return count;
    }

    private static String architectureComponents(
            boolean hasTutr,
            boolean hasAle
    ) {
        return ApaArchitecture.components(hasTutr, hasAle);
    }

    private static void accumulateSample(
            GeneAccumulator accumulator,
            List<SiteCluster> geneSites,
            Sample sample,
            String control,
            String caseCondition
    ) {
        int[] compartmentCounts = new int[COMPARTMENT_COUNT];
        int total = 0;

        for (SiteCluster site : geneSites) {
            int siteCount = count(site, sample.name);
            total += siteCount;
            compartmentCounts[compartment(site)] += siteCount;
        }

        accumulateCompartmentFractions(
                accumulator,
                sample,
                compartmentCounts,
                total,
                control,
                caseCondition
        );
        accumulateBranchFractions(
                accumulator,
                sample,
                compartmentCounts[TERMINAL],
                control,
                caseCondition
        );
        accumulateTutrFractions(
                accumulator,
                sample,
                control,
                caseCondition
        );
    }

    private static void accumulateCompartmentFractions(
            GeneAccumulator accumulator,
            Sample sample,
            int[] counts,
            int total,
            String control,
            String caseCondition
    ) {
        if (total <= 0) {
            return;
        }

        double[] fractions = new double[COMPARTMENT_COUNT];
        for (int i = 0; i < COMPARTMENT_COUNT; i++) {
            fractions[i] = (double) counts[i] / total;
        }

        if (sample.condition.equals(control)) {
            addInPlace(accumulator.controlCompartments, fractions);
            accumulator.controlCompartmentSamples++;
        } else if (sample.condition.equals(caseCondition)) {
            addInPlace(accumulator.caseCompartments, fractions);
            accumulator.caseCompartmentSamples++;
        }
    }

    private static void accumulateBranchFractions(
            GeneAccumulator accumulator,
            Sample sample,
            int terminalCount,
            String control,
            String caseCondition
    ) {
        if (terminalCount <= 0 || accumulator.anchors.isEmpty()) {
            return;
        }

        double[] branchFractions = new double[accumulator.anchors.size()];
        for (int branchIndex = 0;
             branchIndex < accumulator.anchors.size();
             branchIndex++) {
            int branchCount = 0;
            for (SiteCluster site : accumulator.branches.get(
                    accumulator.anchors.get(branchIndex))) {
                branchCount += count(site, sample.name);
            }
            branchFractions[branchIndex] = (double) branchCount / terminalCount;
        }

        if (sample.condition.equals(control)) {
            addInPlace(accumulator.controlBranches, branchFractions);
            accumulator.controlBranchSamples++;
        } else if (sample.condition.equals(caseCondition)) {
            addInPlace(accumulator.caseBranches, branchFractions);
            accumulator.caseBranchSamples++;
        }
    }

    private static void accumulateTutrFractions(
            GeneAccumulator accumulator,
            Sample sample,
            String control,
            String caseCondition
    ) {
        for (int branchIndex = 0;
             branchIndex < accumulator.anchors.size();
             branchIndex++) {
            List<SiteCluster> branch = accumulator.branches.get(
                    accumulator.anchors.get(branchIndex)
            );
            if (branch.size() < 2) {
                continue;
            }

            int branchTotal = 0;
            for (SiteCluster site : branch) {
                branchTotal += count(site, sample.name);
            }
            if (branchTotal == 0) {
                continue;
            }

            List<SiteCluster> orderedSites = new ArrayList<SiteCluster>(branch);
            sortInTranscriptOrder(orderedSites, accumulator.strand);
            double[] samplePau = new double[orderedSites.size()];
            double sampleWpi = 0.0;

            for (int siteIndex = 0; siteIndex < orderedSites.size(); siteIndex++) {
                double normalizedPosition =
                        (double) siteIndex / (orderedSites.size() - 1);
                double pau = (double) count(
                        orderedSites.get(siteIndex),
                        sample.name
                ) / branchTotal;
                samplePau[siteIndex] = pau;
                sampleWpi += normalizedPosition * pau;
            }

            if (sample.condition.equals(control)) {
                accumulator.controlTutrWpi[branchIndex] += sampleWpi;
                accumulator.controlTutrSamples[branchIndex]++;
                addInPlace(accumulator.controlTutrPau.get(branchIndex), samplePau);
            } else if (sample.condition.equals(caseCondition)) {
                accumulator.caseTutrWpi[branchIndex] += sampleWpi;
                accumulator.caseTutrSamples[branchIndex]++;
                addInPlace(accumulator.caseTutrPau.get(branchIndex), samplePau);
            }
        }
    }

    private static void finalizeMetrics(Metrics metrics, GeneAccumulator accumulator) {
        mean(accumulator.controlCompartments, accumulator.controlCompartmentSamples);
        mean(accumulator.caseCompartments, accumulator.caseCompartmentSamples);
        mean(accumulator.controlBranches, accumulator.controlBranchSamples);
        mean(accumulator.caseBranches, accumulator.caseBranchSamples);

        fillCompartmentMetrics(metrics, accumulator);
        fillAleMetrics(metrics, accumulator);
        fillTutrMetrics(metrics, accumulator);
    }

    private static void fillCompartmentMetrics(
            Metrics metrics,
            GeneAccumulator accumulator
    ) {
        double[] control = accumulator.controlCompartments;
        double[] caseValues = accumulator.caseCompartments;

        metrics.controlTerminalUsage = control[TERMINAL];
        metrics.caseTerminalUsage = caseValues[TERMINAL];
        metrics.deltaTerminalUsage = caseValues[TERMINAL] - control[TERMINAL];

        metrics.controlIpaUsage = control[IPA];
        metrics.caseIpaUsage = caseValues[IPA];
        metrics.deltaIpaUsage = caseValues[IPA] - control[IPA];

        metrics.controlInternalUsage = control[INTERNAL];
        metrics.caseInternalUsage = caseValues[INTERNAL];
        metrics.deltaInternalUsage = caseValues[INTERNAL] - control[INTERNAL];

        metrics.controlDownstreamUsage = control[DOWNSTREAM];
        metrics.caseDownstreamUsage = caseValues[DOWNSTREAM];
        metrics.deltaDownstreamUsage = caseValues[DOWNSTREAM] - control[DOWNSTREAM];

        metrics.terminationShiftIndex = totalVariationDistance(control, caseValues);
        metrics.terminationJsd = jensenShannonDivergence(control, caseValues);
        metrics.dominantControl = COMPARTMENTS[argmax(control)];
        metrics.dominantCase = COMPARTMENTS[argmax(caseValues)];
        metrics.transition = metrics.dominantControl + "->" + metrics.dominantCase;
    }

    private static void fillAleMetrics(
            Metrics metrics,
            GeneAccumulator accumulator
    ) {
        if (accumulator.anchors.size() < 2) {
            metrics.aleSwitchIndex = 0.0;
            metrics.aleJsd = 0.0;
            metrics.controlAleWpi = Double.NaN;
            metrics.caseAleWpi = Double.NaN;
            metrics.deltaAleWpi = Double.NaN;
            return;
        }

        metrics.aleSwitchIndex = totalVariationDistance(
                accumulator.controlBranches,
                accumulator.caseBranches
        );
        metrics.aleJsd = jensenShannonDivergence(
                accumulator.controlBranches,
                accumulator.caseBranches
        );
        metrics.controlAleWpi = weightedPositionIndex(accumulator.controlBranches);
        metrics.caseAleWpi = weightedPositionIndex(accumulator.caseBranches);
        metrics.deltaAleWpi = metrics.caseAleWpi - metrics.controlAleWpi;
    }

    private static void fillTutrMetrics(
            Metrics metrics,
            GeneAccumulator accumulator
    ) {
        double totalWeight = 0.0;
        double weightedControlWpi = 0.0;
        double weightedCaseWpi = 0.0;
        double weightedShiftMagnitude = 0.0;

        for (int branchIndex = 0;
             branchIndex < accumulator.anchors.size();
             branchIndex++) {
            List<SiteCluster> branch = accumulator.branches.get(
                    accumulator.anchors.get(branchIndex)
            );
            if (branch.size() < 2
                    || accumulator.controlTutrSamples[branchIndex] == 0
                    || accumulator.caseTutrSamples[branchIndex] == 0) {
                continue;
            }

            double controlWpi = accumulator.controlTutrWpi[branchIndex]
                    / accumulator.controlTutrSamples[branchIndex];
            double caseWpi = accumulator.caseTutrWpi[branchIndex]
                    / accumulator.caseTutrSamples[branchIndex];

            double branchWeight = 0.5 * accumulator.controlBranches[branchIndex]
                    + 0.5 * accumulator.caseBranches[branchIndex];
            if (branchWeight <= 0.0) {
                branchWeight = 1.0;
            }

            double[] controlPau = accumulator.controlTutrPau.get(branchIndex).clone();
            double[] casePau = accumulator.caseTutrPau.get(branchIndex).clone();
            mean(controlPau, accumulator.controlTutrSamples[branchIndex]);
            mean(casePau, accumulator.caseTutrSamples[branchIndex]);

            double branchShiftMagnitude = totalVariationDistance(controlPau, casePau);
            weightedControlWpi += branchWeight * controlWpi;
            weightedCaseWpi += branchWeight * caseWpi;
            weightedShiftMagnitude += branchWeight * branchShiftMagnitude;
            totalWeight += branchWeight;
        }

        if (totalWeight > 0.0) {
            metrics.controlTutrWpi = weightedControlWpi / totalWeight;
            metrics.caseTutrWpi = weightedCaseWpi / totalWeight;
            metrics.deltaTutrWpi = metrics.caseTutrWpi - metrics.controlTutrWpi;
            metrics.tutrShiftMagnitude = weightedShiftMagnitude / totalWeight;
        } else {
            metrics.controlTutrWpi = Double.NaN;
            metrics.caseTutrWpi = Double.NaN;
            metrics.deltaTutrWpi = Double.NaN;
            metrics.tutrShiftMagnitude = 0.0;
        }
    }

    public static void enrichGeneResults(
            List<DifferentialApa.GeneResult> genes,
            Map<String, Metrics> metricsByGene,
            double fdr,
            double internalExonThreshold,
            double aleThreshold,
            double ipaThreshold,
            double tutrThreshold
    ) {
        for (DifferentialApa.GeneResult gene : genes) {
            Metrics metrics = metricsByGene.get(gene.geneId);
            if (metrics == null) {
                continue;
            }

            copyMetrics(gene, metrics);
            classifyGeneChange(
                    gene,
                    metrics,
                    fdr,
                    internalExonThreshold,
                    aleThreshold,
                    ipaThreshold,
                    tutrThreshold
            );
        }
    }

    private static void copyMetrics(
            DifferentialApa.GeneResult gene,
            Metrics metrics
    ) {
        gene.apaType = metrics.architectureClass;
        gene.architectureComponents = metrics.components;
        gene.nTerminalBranches = metrics.nTerminalBranches;
        gene.nTutrBranches = metrics.nTutrBranches;
        gene.controlTerminalUsage = metrics.controlTerminalUsage;
        gene.caseTerminalUsage = metrics.caseTerminalUsage;
        gene.deltaTerminalUsage = metrics.deltaTerminalUsage;
        gene.controlIpaUsage = metrics.controlIpaUsage;
        gene.caseIpaUsage = metrics.caseIpaUsage;
        gene.deltaIpaUsage = metrics.deltaIpaUsage;
        gene.terminationShiftIndex = metrics.terminationShiftIndex;
        gene.terminationJsd = metrics.terminationJsd;
        gene.aleSwitchIndex = metrics.aleSwitchIndex;
        gene.aleJsd = metrics.aleJsd;
        gene.controlAleWpi = metrics.controlAleWpi;
        gene.caseAleWpi = metrics.caseAleWpi;
        gene.deltaAleWpi = metrics.deltaAleWpi;
        gene.controlTutrWpi = metrics.controlTutrWpi;
        gene.caseTutrWpi = metrics.caseTutrWpi;
        gene.deltaTutrWpi = metrics.deltaTutrWpi;
        gene.tutrShiftMagnitude = metrics.tutrShiftMagnitude;
        gene.dominantControlArchitecture = metrics.dominantControl;
        gene.dominantCaseArchitecture = metrics.dominantCase;
        gene.architectureTransition = metrics.transition;
    }

    private static void classifyGeneChange(
            DifferentialApa.GeneResult gene,
            Metrics metrics,
            double fdr,
            double internalExonThreshold,
            double aleThreshold,
            double ipaThreshold,
            double tutrThreshold
    ) {
        boolean supported = gene.ksFdr <= fdr || gene.minFdr <= fdr;
        initializeDirections(gene, metrics);

        List<String> activeDirections = new ArrayList<String>();
        List<String> changedMechanisms = new ArrayList<String>();
        List<String> mechanismDirections = new ArrayList<String>();

        if (supported) {
            classifyIpa(
                    gene,
                    metrics,
                    fdr,
                    ipaThreshold,
                    activeDirections,
                    changedMechanisms,
                    mechanismDirections
            );
            classifyAle(
                    gene,
                    metrics,
                    fdr,
                    aleThreshold,
                    tutrThreshold,
                    activeDirections,
                    changedMechanisms,
                    mechanismDirections
            );
            classifyTutr(
                    gene,
                    metrics,
                    fdr,
                    tutrThreshold,
                    activeDirections,
                    changedMechanisms,
                    mechanismDirections
            );
            classifyInternalExon(
                    gene,
                    metrics,
                    fdr,
                    internalExonThreshold,
                    activeDirections,
                    changedMechanisms,
                    mechanismDirections
            );
        }

        gene.apaChange = supported && !changedMechanisms.isEmpty();
        gene.apaDirection = overallDirection(supported, gene.apaChange, activeDirections);
        gene.changedMechanisms = changedMechanisms.isEmpty()
                ? "none"
                : join(changedMechanisms, "+");
        gene.mechanismDirections = mechanismDirections.isEmpty()
                ? "none"
                : join(mechanismDirections, ";");
        gene.apaEvent = buildApaEvent(gene, changedMechanisms);
    }

    private static void initializeDirections(
            DifferentialApa.GeneResult gene,
            Metrics metrics
    ) {
        gene.ipaDirection = metrics.controlIpaUsage > 0 || metrics.caseIpaUsage > 0
                ? "stable"
                : "not_applicable";
        gene.aleDirection = metrics.nTerminalBranches >= 2
                ? "stable"
                : "not_applicable";
        gene.tutrDirection = metrics.nTutrBranches > 0
                ? "stable"
                : "not_applicable";
        gene.internalExonDirection =
                metrics.controlInternalUsage > 0 || metrics.caseInternalUsage > 0
                        ? "stable"
                        : "not_applicable";
    }

    private static void classifyIpa(
            DifferentialApa.GeneResult gene,
            Metrics metrics,
            double fdr,
            double threshold,
            List<String> activeDirections,
            List<String> changedMechanisms,
            List<String> mechanismDirections
    ) {
        if (gene.ipaMinFdr > fdr
                || "not_applicable".equals(gene.ipaDirection)
                || Math.abs(metrics.deltaIpaUsage) < threshold) {
            return;
        }

        gene.ipaDirection = metrics.deltaIpaUsage > 0 ? "gain" : "loss";
        addMechanism(
                "IPA",
                gene.ipaDirection,
                activeDirections,
                changedMechanisms,
                mechanismDirections
        );
    }

    private static void classifyAle(
            DifferentialApa.GeneResult gene,
            Metrics metrics,
            double fdr,
            double switchThreshold,
            double directionThreshold,
            List<String> activeDirections,
            List<String> changedMechanisms,
            List<String> mechanismDirections
    ) {
        if (gene.aleMinFdr > fdr
                || "not_applicable".equals(gene.aleDirection)
                || metrics.aleSwitchIndex < switchThreshold) {
            return;
        }

        if (!Double.isNaN(metrics.deltaAleWpi)
                && Math.abs(metrics.deltaAleWpi) >= directionThreshold) {
            gene.aleDirection = metrics.deltaAleWpi > 0
                    ? "distal_switch"
                    : "proximal_switch";
        } else {
            gene.aleDirection = "switch";
        }

        addMechanism(
                "ALE",
                gene.aleDirection,
                activeDirections,
                changedMechanisms,
                mechanismDirections
        );
    }

    private static void classifyTutr(
            DifferentialApa.GeneResult gene,
            Metrics metrics,
            double fdr,
            double threshold,
            List<String> activeDirections,
            List<String> changedMechanisms,
            List<String> mechanismDirections
    ) {
        if (gene.tutrMinFdr > fdr
                || "not_applicable".equals(gene.tutrDirection)
                || metrics.tutrShiftMagnitude < threshold) {
            return;
        }

        if (!Double.isNaN(metrics.deltaTutrWpi)
                && Math.abs(metrics.deltaTutrWpi) >= threshold) {
            gene.tutrDirection = metrics.deltaTutrWpi > 0
                    ? "lengthening"
                    : "shortening";
        } else {
            gene.tutrDirection = "no_directional_shift";
        }

        addMechanism(
                "TUTR",
                gene.tutrDirection,
                activeDirections,
                changedMechanisms,
                mechanismDirections
        );
    }

    private static void classifyInternalExon(
            DifferentialApa.GeneResult gene,
            Metrics metrics,
            double fdr,
            double threshold,
            List<String> activeDirections,
            List<String> changedMechanisms,
            List<String> mechanismDirections
    ) {
        if (gene.internalExonMinFdr > fdr
                || "not_applicable".equals(gene.internalExonDirection)
                || Math.abs(metrics.deltaInternalUsage) < threshold) {
            return;
        }

        gene.internalExonDirection = metrics.deltaInternalUsage > 0
                ? "gain"
                : "loss";
        addMechanism(
                "internal_exon",
                gene.internalExonDirection,
                activeDirections,
                changedMechanisms,
                mechanismDirections
        );
    }

    private static void addMechanism(
            String mechanism,
            String direction,
            List<String> activeDirections,
            List<String> changedMechanisms,
            List<String> mechanismDirections
    ) {
        activeDirections.add(direction);
        changedMechanisms.add(mechanism);
        mechanismDirections.add(mechanism + "=" + direction);
    }

    private static String overallDirection(
            boolean supported,
            boolean apaChange,
            List<String> activeDirections
    ) {
        if (!supported) {
            return "stable";
        }
        if (activeDirections.size() > 1) {
            return "complex";
        }
        if (activeDirections.size() == 1) {
            return activeDirections.get(0);
        }
        if (apaChange) {
            return "no_directional_shift";
        }
        return "stable";
    }

    private static String buildApaEvent(
            DifferentialApa.GeneResult gene,
            List<String> changedMechanisms
    ) {
        if (changedMechanisms.size() == 1) {
            return changedMechanisms.get(0) + ":" + gene.apaDirection;
        }
        if (changedMechanisms.size() > 1) {
            return "mixed".equals(gene.apaType)
                    ? "mixed:complex"
                    : "multi_mechanism:complex";
        }
        return gene.apaType + ":" + gene.apaDirection;
    }

    private static boolean hasUnresolvedTerminalPas(
            List<SiteCluster> geneSites
    ) {
        for (SiteCluster site : geneSites) {
            if (site.terminalExon && site.terminalStructureAnchor0 < 0) {
                return true;
            }
        }
        return false;
    }

    private static String architectureClass(
            int terminalStructures,
            int tutrStructures,
            int terminalPas
    ) {
        return ApaArchitecture.classify(
                terminalStructures,
                tutrStructures,
                terminalPas
        );
    }

    private static Map<String, List<SiteCluster>> groupByGene(List<SiteCluster> sites) {
        Map<String, List<SiteCluster>> grouped =
                new LinkedHashMap<String, List<SiteCluster>>();

        for (SiteCluster site : sites) {
            if (".".equals(site.geneId)) {
                continue;
            }
            List<SiteCluster> geneSites = grouped.get(site.geneId);
            if (geneSites == null) {
                geneSites = new ArrayList<SiteCluster>();
                grouped.put(site.geneId, geneSites);
            }
            geneSites.add(site);
        }
        return grouped;
    }

    private static Map<Integer, List<SiteCluster>> terminalBranches(
            List<SiteCluster> geneSites
    ) {
        Map<Integer, List<SiteCluster>> branches =
                new LinkedHashMap<Integer, List<SiteCluster>>();

        for (SiteCluster site : geneSites) {
            if (!site.terminalExon
                    || site.terminalStructureAnchor0 < 0) {
                continue;
            }

            List<SiteCluster> branchSites =
                    branches.get(
                            site.terminalStructureAnchor0
                    );

            if (branchSites == null) {
                branchSites =
                        new ArrayList<SiteCluster>();

                branches.put(
                        site.terminalStructureAnchor0,
                        branchSites
                );
            }

            branchSites.add(site);
        }

        return branches;
    }

    private static List<Integer> orderedAnchors(
            final Map<Integer, List<SiteCluster>> branches,
            final Strand strand
    ) {
        List<Integer> anchors = new ArrayList<Integer>(branches.keySet());
        Collections.sort(anchors, new Comparator<Integer>() {
            @Override
            public int compare(Integer left, Integer right) {
                int comparison = Double.compare(
                        branchCenter(branches.get(left)),
                        branchCenter(branches.get(right))
                );
                return strand == Strand.MINUS ? -comparison : comparison;
            }
        });
        return anchors;
    }

    private static double branchCenter(List<SiteCluster> sites) {
        double total = 0.0;
        for (SiteCluster site : sites) {
            total += site.representativePos0;
        }
        return sites.isEmpty() ? 0.0 : total / sites.size();
    }

    private static void sortInTranscriptOrder(
            List<SiteCluster> sites,
            final Strand strand
    ) {
        Collections.sort(sites, new Comparator<SiteCluster>() {
            @Override
            public int compare(SiteCluster left, SiteCluster right) {
                int comparison = Integer.compare(
                        left.representativePos0,
                        right.representativePos0
                );
                return strand == Strand.MINUS ? -comparison : comparison;
            }
        });
    }

    private static int compartment(SiteCluster site) {
        if (site.terminalExon) {
            return TERMINAL;
        }
        if ("INTRON".equals(site.region)) {
            return IPA;
        }
        if ("EXON_CODING".equals(site.region)
                || "EXON_NONCODING".equals(site.region)) {
            return INTERNAL;
        }
        return DOWNSTREAM;
    }

    private static int count(SiteCluster site, String sampleName) {
        Integer value = site.counts.get(sampleName);
        return value == null ? 0 : value;
    }

    private static void addInPlace(double[] target, double[] values) {
        for (int i = 0; i < Math.min(target.length, values.length); i++) {
            target[i] += values[i];
        }
    }

    private static void mean(double[] values, int count) {
        if (count <= 0) {
            return;
        }
        for (int i = 0; i < values.length; i++) {
            values[i] /= count;
        }
    }

    private static double totalVariationDistance(double[] left, double[] right) {
        double distance = 0.0;
        for (int i = 0; i < Math.min(left.length, right.length); i++) {
            distance += Math.abs(left[i] - right[i]);
        }
        return 0.5 * distance;
    }

    private static double jensenShannonDivergence(double[] left, double[] right) {
        if (left.length == 0 || right.length == 0) {
            return 0.0;
        }

        double leftTotal = sum(left);
        double rightTotal = sum(right);
        if (leftTotal <= 0.0 || rightTotal <= 0.0) {
            return 0.0;
        }

        double divergence = 0.0;
        for (int i = 0; i < Math.min(left.length, right.length); i++) {
            double p = left[i] / leftTotal;
            double q = right[i] / rightTotal;
            double midpoint = 0.5 * (p + q);

            if (p > 0.0) {
                divergence += 0.5 * p * log2(p / midpoint);
            }
            if (q > 0.0) {
                divergence += 0.5 * q * log2(q / midpoint);
            }
        }
        return Math.max(0.0, Math.min(1.0, divergence));
    }

    private static double log2(double value) {
        return Math.log(value) / Math.log(2.0);
    }

    private static double sum(double[] values) {
        double total = 0.0;
        for (double value : values) {
            total += value;
        }
        return total;
    }

    private static int argmax(double[] values) {
        int bestIndex = 0;
        for (int i = 1; i < values.length; i++) {
            if (values[i] > values[bestIndex]) {
                bestIndex = i;
            }
        }
        return bestIndex;
    }

    private static double weightedPositionIndex(double[] probabilities) {
        if (probabilities.length == 0) {
            return Double.NaN;
        }

        double total = sum(probabilities);
        if (total <= 0.0) {
            return Double.NaN;
        }
        if (probabilities.length == 1) {
            return 0.5;
        }

        double weightedPosition = 0.0;
        for (int i = 0; i < probabilities.length; i++) {
            double normalizedPosition = (double) i / (probabilities.length - 1);
            weightedPosition += (probabilities[i] / total) * normalizedPosition;
        }
        return weightedPosition;
    }

    private static String join(List<String> values, String separator) {
        if (values.isEmpty()) {
            return "none";
        }

        StringBuilder joined = new StringBuilder();
        for (String value : values) {
            if (joined.length() > 0) {
                joined.append(separator);
            }
            joined.append(value);
        }
        return joined.toString();
    }
}
