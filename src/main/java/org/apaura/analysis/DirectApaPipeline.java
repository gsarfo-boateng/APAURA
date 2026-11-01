package org.apaura.analysis;

import org.apaura.annotation.GtfAnnotation;
import org.apaura.io.SampleSheet;
import org.apaura.model.EndKey;
import org.apaura.model.Sample;
import org.apaura.model.SiteCluster;
import org.apaura.output.CanonicalOutputWriter;
import org.apaura.output.LongReadCanonicalAdapter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/** End-to-end long-read/direct-RNA APA analysis pipeline. */
public final class DirectApaPipeline {
    private DirectApaPipeline() {
    }

    public static void run(
            String sampleSheet,
            String gtf,
            String outDir,
            String control,
            String caseCondition,
            int minMapq,
            boolean keepDuplicates,
            String strandMode,
            int clusterDistance,
            int annotationSlop,
            int minTotal,
            double minDeltaPau,
            double fdr,
            double wpiThreshold
    ) throws Exception {
        run(
                sampleSheet,
                gtf,
                outDir,
                control,
                caseCondition,
                minMapq,
                keepDuplicates,
                strandMode,
                clusterDistance,
                annotationSlop,
                minTotal,
                minDeltaPau,
                fdr,
                wpiThreshold,
                "bam",
                "mode",
                10,
                5,
                0.01,
                0.1,
                0.1,
                0.1,
                0.1
        );
    }

    public static void run(
            String sampleSheet,
            String gtf,
            String outDir,
            String control,
            String caseCondition,
            int minMapq,
            boolean keepDuplicates,
            String strandMode,
            int clusterDistance,
            int annotationSlop,
            int minTotal,
            double minDeltaPau,
            double fdr,
            double wpiThreshold,
            String inputFormat,
            String clusterRepresentative,
            int minGeneReads,
            int minSiteReads,
            double minPau,
            double internalExonThreshold,
            double aleThreshold,
            double ipaThreshold,
            double tutrThreshold
    ) throws Exception {
        File outputDirectory = ensureOutputDirectory(outDir);
        File detailsDirectory =
                CanonicalOutputWriter.prepareDetailsDirectory(
                        outputDirectory
                );
        List<Sample> samples = SampleSheet.read(sampleSheet);
        validateConditions(samples, control, caseCondition);

        Map<String, Map<EndKey, Integer>> endCounts = extractEnds(
                samples,
                inputFormat,
                minMapq,
                keepDuplicates,
                strandMode
        );

        List<SiteCluster> allSites = SiteClusterer.cluster(
                samples,
                endCounts,
                clusterDistance,
                clusterRepresentative
        );
        System.err.println(
                "[APAURA] pooled PAS clusters="
                        + allSites.size()
                        + " cluster_distance="
                        + clusterDistance
                        + " representative="
                        + clusterRepresentative
        );

        if (gtf == null) {
            throw new IllegalArgumentException(
                    "Differential long-read analysis requires --gtf for gene/PAS annotation"
            );
        }

        GtfAnnotation annotation = GtfAnnotation.load(gtf);
        annotation.annotate(allSites, annotationSlop);

        TerminalStructureResolver.assign(
                annotation.transcripts(),
                allSites
        );

        writeSites(
                new File(detailsDirectory, "pas_all.bed"),
                allSites
        );

        DirectSiteFilter.Result filtered = DirectSiteFilter.filter(
                samples,
                allSites,
                control,
                caseCondition,
                minGeneReads,
                minSiteReads,
                minPau
        );
        List<SiteCluster> retainedSites = filtered.retained;
        ArchitectureAnalysis.annotateSiteRoles(retainedSites);

        writeSites(new File(detailsDirectory, "pas.bed"), retainedSites);
        writeMatrix(new File(detailsDirectory, "pas_counts.tsv"), samples, retainedSites);
        writeUsage(
                new File(detailsDirectory, "pau.tsv"),
                samples,
                retainedSites,
                filtered.geneTotals
        );

        DifferentialApa.Results differential = DifferentialApa.run(
                samples,
                retainedSites,
                control,
                caseCondition,
                minTotal,
                minDeltaPau,
                fdr,
                wpiThreshold,
                filtered.geneTotals
        );

        Map<String, ArchitectureAnalysis.Metrics> architecture =
                ArchitectureAnalysis.analyze(
                        samples,
                        retainedSites,
                        control,
                        caseCondition
                );
        ArchitectureAnalysis.enrichGeneResults(
                differential.genes,
                architecture,
                fdr,
                internalExonThreshold,
                aleThreshold,
                ipaThreshold,
                tutrThreshold
        );

        
        // APAURA_RUN_3SEQ_DM_GATE
        if (inputFormat != null
                && inputFormat.toLowerCase(Locale.ROOT).startsWith("3seq:")) {

            ThreeSeqGeneDifferential.Batch threeSeqDm =
                    ThreeSeqGeneDifferential.test(
                            samples,
                            retainedSites,
                            control,
                            caseCondition
                    );

            ThreeSeqGeneDifferential.applyGate(
                    differential.genes,
                    threeSeqDm,
                    fdr
            );

            File threeSeqDmFile =
                    new File(
                            new File(outputDirectory, "details"),
                            "three_seq_dm.tsv"
                    );

            ThreeSeqGeneDifferential.write(
                    threeSeqDmFile,
                    threeSeqDm
            );

            System.err.println(
                    "[APAURA] run-3seq gene-level moderated DM: "
                            + "common_precision="
                            + threeSeqDm.commonPrecision
                            + " output="
                            + threeSeqDmFile.getAbsolutePath()
            );
        }

writeSiteDiff(new File(detailsDirectory, "apa_sites.tsv"), differential.sites);
        writeGene(new File(detailsDirectory, "apa.tsv"), differential.genes);
        writeArchitecture(
                new File(detailsDirectory, "apa_architecture.tsv"),
                differential.genes
        );
        writeRunInfo(
                new File(detailsDirectory, "run_info.txt"),
                inputFormat,
                clusterDistance,
                clusterRepresentative,
                minGeneReads,
                minSiteReads,
                minPau,
                minDeltaPau,
                fdr,
                wpiThreshold,
                internalExonThreshold,
                aleThreshold,
                ipaThreshold,
                tutrThreshold,
                control,
                caseCondition
        );

        LongReadCanonicalAdapter.write(
                outputDirectory,
                sampleSheet,
                gtf,
                inputFormat,
                control,
                caseCondition,
                clusterDistance,
                clusterRepresentative,
                minGeneReads,
                minSiteReads,
                minPau,
                minDeltaPau,
                fdr,
                wpiThreshold,
                internalExonThreshold,
                aleThreshold,
                ipaThreshold,
                tutrThreshold,
                retainedSites,
                differential,
                architecture
        );
    }

    private static File ensureOutputDirectory(String outDir) throws IOException {
        File directory = new File(outDir);
        if (!directory.exists() && !directory.mkdirs()) {
            throw new IOException("Cannot create output directory " + outDir);
        }
        return directory;
    }

    private static void validateConditions(
            List<Sample> samples,
            String control,
            String caseCondition
    ) {
        Set<String> conditions = new LinkedHashSet<String>();
        for (Sample sample : samples) {
            conditions.add(sample.condition);
        }

        if (!conditions.contains(control) || !conditions.contains(caseCondition)) {
            throw new IllegalArgumentException(
                    "Conditions in sample sheet="
                            + conditions
                            + "; requested control="
                            + control
                            + " case="
                            + caseCondition
            );
        }
    }

    private static Map<String, Map<EndKey, Integer>> extractEnds(
            List<Sample> samples,
            String inputFormat,
            int minMapq,
            boolean keepDuplicates,
            String strandMode
    ) throws Exception {
        Map<String, Map<EndKey, Integer>> endCounts =
                new LinkedHashMap<String, Map<EndKey, Integer>>();

        for (Sample sample : samples) {
            Map<EndKey, Integer> sampleEnds;
            if ("isoquant".equalsIgnoreCase(inputFormat)) {
                sampleEnds = IsoQuantEndExtractor.extract(sample.path, strandMode);
            } else if (
                    inputFormat != null
                            && inputFormat.toLowerCase(Locale.ROOT).startsWith("3seq:")
            ) {
                // APAURA_RUN_3SEQ_EXTRACT
                String protocolName = inputFormat.substring("3seq:".length());
                ProtocolSpec protocol = ProtocolRegistry.get(protocolName);

                sampleEnds = EndExtractor.extract(
                        sample.path,
                        minMapq,
                        keepDuplicates,
                        strandMode,
                        protocol
                );

            } else {
                sampleEnds = EndExtractor.extract(
                        sample.path,
                        minMapq,
                        keepDuplicates,
                        strandMode
                );
            }
            endCounts.put(sample.name, sampleEnds);
        }
        return endCounts;
    }

    private static void writeSites(File file, List<SiteCluster> sites) throws IOException {
        try (PrintWriter writer = new PrintWriter(
                new BufferedWriter(new FileWriter(file)))) {
            writer.println(
                    "#chrom\tstart\tend\tsite_id\ttotal_count\tstrand\tcluster_start\t"
                            + "cluster_end\tgene_id\tgene_name\tbiotype\tregion\t"
                            + "architecture_role\tterminal_exon\tterminal_exon_anchor0\t"
                            + "terminal_branch_rank\tsite_rank\tnormalized_position"
            );

            for (SiteCluster site : sites) {
                writer.printf(
                        Locale.ROOT,
                        "%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s%n",
                        site.chrom,
                        site.representativePos0,
                        site.representativePos0 + 1,
                        site.id,
                        site.totalCount(),
                        site.strand.symbol(),
                        site.startPos0,
                        site.endPos0 + 1,
                        site.geneId,
                        site.geneName,
                        site.biotype,
                        displayRegion(site.region),
                        site.architectureRole,
                        site.terminalExon ? "true" : "false",
                        optionalInt(site.terminalExonAnchor0),
                        optionalInt(site.terminalBranchRank),
                        site.siteRank,
                        Double.isNaN(site.normalizedPosition)
                                ? "NA"
                                : String.format(Locale.ROOT, "%.6f", site.normalizedPosition)
                );
            }
        }
    }

    private static void writeMatrix(
            File file,
            List<Sample> samples,
            List<SiteCluster> sites
    ) throws IOException {
        try (PrintWriter writer = new PrintWriter(
                new BufferedWriter(new FileWriter(file)))) {
            writer.print("site_id\tchrom\tpos0\tstrand\tgene_id\tgene_name\tregion");
            for (Sample sample : samples) {
                writer.print("\t" + sample.name);
            }
            writer.println();

            for (SiteCluster site : sites) {
                writer.printf(
                        "%s\t%s\t%d\t%c\t%s\t%s\t%s",
                        site.id,
                        site.chrom,
                        site.representativePos0,
                        site.strand.symbol(),
                        site.geneId,
                        site.geneName,
                        displayRegion(site.region)
                );
                for (Sample sample : samples) {
                    writer.print("\t" + count(site, sample.name));
                }
                writer.println();
            }
        }
    }

    private static void writeUsage(
            File file,
            List<Sample> samples,
            List<SiteCluster> sites,
            Map<String, Map<String, Integer>> geneTotals
    ) throws IOException {
        try (PrintWriter writer = new PrintWriter(
                new BufferedWriter(new FileWriter(file)))) {
            writer.println(
                    "sample\tcondition\tsite_id\tgene_id\tgene_name\tchrom\tpos0\tstrand\t"
                            + "region\tsite_rank\tcount\tgene_end_count\tpau"
            );

            for (SiteCluster site : sites) {
                for (Sample sample : samples) {
                    int siteCount = count(site, sample.name);
                    int geneTotal = geneTotal(geneTotals, site.geneId, sample.name);
                    double pau = geneTotal > 0
                            ? (double) siteCount / geneTotal
                            : Double.NaN;

                    writer.printf(
                            Locale.ROOT,
                            "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%c\t%s\t%d\t%d\t%d\t%s%n",
                            sample.name,
                            sample.condition,
                            site.id,
                            site.geneId,
                            site.geneName,
                            site.chrom,
                            site.representativePos0,
                            site.strand.symbol(),
                            displayRegion(site.region),
                            site.siteRank,
                            siteCount,
                            geneTotal,
                            Double.isNaN(pau)
                                    ? "NA"
                                    : String.format(Locale.ROOT, "%.8f", pau)
                    );
                }
            }
        }
    }

    private static void writeSiteDiff(
            File file,
            List<DifferentialApa.SiteResult> results
    ) throws IOException {
        try (PrintWriter writer = new PrintWriter(
                new BufferedWriter(new FileWriter(file)))) {
            writer.println(
                    "site_id\tgene_id\tgene_name\tbiotype\tchrom\tpos0\tstrand\tregion\t"
                            + "architecture_role\tterminal_exon\tterminal_exon_anchor0\t"
                            + "terminal_branch_rank\tsite_rank\tcontrol_pau\tcase_pau\t"
                            + "delta_pau\tcontrol_n_evaluable\tcase_n_evaluable\t"
                            + "pvalue\tfdr\tdispersion_rho\tmethod\tstatus"
            );

            for (DifferentialApa.SiteResult result : results) {
                SiteCluster site = result.site;
                writer.printf(
                        Locale.ROOT,
                        "%s\t%s\t%s\t%s\t%s\t%d\t%c\t%s\t%s\t%s\t%s\t%s\t%d\t%.6f\t%.6f\t%.6f\t%d\t%d\t%s\t%s\t%s\t%s\t%s%n",
                        site.id,
                        site.geneId,
                        site.geneName,
                        site.biotype,
                        site.chrom,
                        site.representativePos0,
                        site.strand.symbol(),
                        displayRegion(site.region),
                        site.architectureRole,
                        site.terminalExon ? "true" : "false",
                        optionalInt(site.terminalExonAnchor0),
                        optionalInt(site.terminalBranchRank),
                        site.siteRank,
                        result.controlPau,
                        result.casePau,
                        result.delta,
                        result.controlReplicates,
                        result.caseReplicates,
                        fmt8(result.p),
                        fmt8(result.q),
                        Double.isNaN(result.rho)
                                ? "NA"
                                : String.format(Locale.ROOT, "%.6g", result.rho),
                        result.method,
                        result.status
                );
            }
        }
    }

    private static void writeGene(
            File file,
            List<DifferentialApa.GeneResult> results
    ) throws IOException {
        try (PrintWriter writer = new PrintWriter(
                new BufferedWriter(new FileWriter(file)))) {
            writer.println(
                    "gene_id\tgene_name\tbiotype\tapa_architecture\tarchitecture_components\t"
                            + "apa_change\tapa_direction\ttutr_direction\tale_direction\t"
                            + "ipa_direction\tinternal_exon_direction\tapa_event\tchanged_mechanisms\t"
                            + "mechanism_directions\tn_sites\tn_terminal_branches\tn_tutr_branches\t"
                            + "control_wpi\tcase_wpi\tdelta_wpi\tapd\tusage_jsd\tcontrol_terminal_usage\t"
                            + "case_terminal_usage\tdelta_terminal_usage\tcontrol_ipa_usage\t"
                            + "case_ipa_usage\tdelta_ipa_usage\ttermination_shift_index\t"
                            + "termination_jsd\tale_switch_index\tale_jsd\tcontrol_ale_wpi\t"
                            + "case_ale_wpi\tdelta_ale_wpi\tcontrol_tutr_wpi\tcase_tutr_wpi\t"
                            + "delta_tutr_wpi\ttutr_shift_magnitude\tdominant_control_termination\t"
                            + "dominant_case_termination\ttermination_transition\tks_d\tks_signed_d\t"
                            + "ks_max_position_1based\tks_z\tks_pvalue\tks_fdr\tmin_site_fdr"
            );

            for (DifferentialApa.GeneResult result : results) {
                writer.println(geneResultLine(result));
            }
        }
    }

    private static String geneResultLine(DifferentialApa.GeneResult result) {
        List<String> fields = new ArrayList<String>();
        fields.add(result.geneId);
        fields.add(result.geneName);
        fields.add(result.biotype);
        fields.add(result.apaType);
        fields.add(result.architectureComponents);
        fields.add(result.apaChange ? "yes" : "no");
        fields.add(result.apaDirection);
        fields.add(result.tutrDirection);
        fields.add(result.aleDirection);
        fields.add(result.ipaDirection);
        fields.add(result.internalExonDirection);
        fields.add(result.apaEvent);
        fields.add(result.changedMechanisms);
        fields.add(result.mechanismDirections);
        fields.add(Integer.toString(result.nSites));
        fields.add(Integer.toString(result.nTerminalBranches));
        fields.add(Integer.toString(result.nTutrBranches));
        fields.add(fmt(result.controlWpi));
        fields.add(fmt(result.caseWpi));
        fields.add(fmt(result.deltaWpi));
        fields.add(fmt(result.apd));
        fields.add(fmt(result.usageJsd));
        fields.add(fmt(result.controlTerminalUsage));
        fields.add(fmt(result.caseTerminalUsage));
        fields.add(fmt(result.deltaTerminalUsage));
        fields.add(fmt(result.controlIpaUsage));
        fields.add(fmt(result.caseIpaUsage));
        fields.add(fmt(result.deltaIpaUsage));
        fields.add(fmt(result.terminationShiftIndex));
        fields.add(fmt(result.terminationJsd));
        fields.add(fmt(result.aleSwitchIndex));
        fields.add(fmt(result.aleJsd));
        fields.add(fmt(result.controlAleWpi));
        fields.add(fmt(result.caseAleWpi));
        fields.add(fmt(result.deltaAleWpi));
        fields.add(fmt(result.controlTutrWpi));
        fields.add(fmt(result.caseTutrWpi));
        fields.add(fmt(result.deltaTutrWpi));
        fields.add(fmt(result.tutrShiftMagnitude));
        fields.add(result.dominantControlArchitecture);
        fields.add(result.dominantCaseArchitecture);
        fields.add(result.architectureTransition);
        fields.add(fmt(result.ksD));
        fields.add(fmt(result.ksSignedD));
        fields.add(result.ksMaxDiffAt0 < 0
                ? "NA"
                : Integer.toString(result.ksMaxDiffAt0 + 1));
        fields.add(fmt(result.ksZ));
        fields.add(fmt8(result.ksPvalue));
        fields.add(fmt8(result.ksFdr));
        fields.add(fmt8(result.minFdr));
        return joinTabs(fields);
    }

    private static void writeArchitecture(
            File file,
            List<DifferentialApa.GeneResult> results
    ) throws IOException {
        try (PrintWriter writer = new PrintWriter(
                new BufferedWriter(new FileWriter(file)))) {
            writer.println(
                    "gene_id\tgene_name\tapa_architecture\tarchitecture_components\tapa_change\t"
                            + "apa_direction\ttutr_direction\tale_direction\tipa_direction\t"
                            + "internal_exon_direction\tapa_event\tchanged_mechanisms\t"
                            + "mechanism_directions\tn_terminal_branches\tn_tutr_branches\t"
                            + "control_terminal_usage\tcase_terminal_usage\tdelta_terminal_usage\t"
                            + "control_ipa_usage\tcase_ipa_usage\tdelta_ipa_usage\t"
                            + "termination_shift_index\ttermination_jsd\tale_switch_index\t"
                            + "delta_ale_wpi\tdelta_tutr_wpi\ttutr_shift_magnitude\t"
                            + "termination_transition"
            );

            for (DifferentialApa.GeneResult result : results) {
                writer.println(architectureResultLine(result));
            }
        }
    }

    private static String architectureResultLine(DifferentialApa.GeneResult result) {
        List<String> fields = new ArrayList<String>();
        fields.add(result.geneId);
        fields.add(result.geneName);
        fields.add(result.apaType);
        fields.add(result.architectureComponents);
        fields.add(result.apaChange ? "yes" : "no");
        fields.add(result.apaDirection);
        fields.add(result.tutrDirection);
        fields.add(result.aleDirection);
        fields.add(result.ipaDirection);
        fields.add(result.internalExonDirection);
        fields.add(result.apaEvent);
        fields.add(result.changedMechanisms);
        fields.add(result.mechanismDirections);
        fields.add(Integer.toString(result.nTerminalBranches));
        fields.add(Integer.toString(result.nTutrBranches));
        fields.add(fmt(result.controlTerminalUsage));
        fields.add(fmt(result.caseTerminalUsage));
        fields.add(fmt(result.deltaTerminalUsage));
        fields.add(fmt(result.controlIpaUsage));
        fields.add(fmt(result.caseIpaUsage));
        fields.add(fmt(result.deltaIpaUsage));
        fields.add(fmt(result.terminationShiftIndex));
        fields.add(fmt(result.terminationJsd));
        fields.add(fmt(result.aleSwitchIndex));
        fields.add(fmt(result.deltaAleWpi));
        fields.add(fmt(result.deltaTutrWpi));
        fields.add(fmt(result.tutrShiftMagnitude));
        fields.add(result.architectureTransition);
        return joinTabs(fields);
    }

    private static int geneTotal(
            Map<String, Map<String, Integer>> geneTotals,
            String geneId,
            String sampleName
    ) {
        Map<String, Integer> sampleTotals = geneTotals.get(geneId);
        if (sampleTotals == null) {
            return 0;
        }
        Integer total = sampleTotals.get(sampleName);
        return total == null ? 0 : total;
    }

    private static String optionalInt(int value) {
        return value < 0 ? "NA" : Integer.toString(value);
    }

    private static String displayRegion(String region) {
        if (region == null) {
            return "unknown";
        }

        switch (region) {
            case "THREE_PRIME_UTR":
                return "3UTR";
            case "TERMINAL_EXON_NONCODING":
                return "terminal_exon_noncoding";
            case "EXON_CODING":
                return "coding_exon";
            case "EXON_NONCODING":
                return "noncoding_exon";
            case "INTRON":
                return "intron";
            case "DOWNSTREAM_3P":
                return "downstream";
            case "INTERGENIC":
                return "intergenic";
            default:
                return region.toLowerCase(Locale.ROOT);
        }
    }

    private static String fmt(double value) {
        return Double.isNaN(value)
                ? "NA"
                : String.format(Locale.ROOT, "%.6f", value);
    }

    private static String fmt8(double value) {
        return Double.isNaN(value)
                ? "NA"
                : String.format(Locale.ROOT, "%.8g", value);
    }

    private static String joinTabs(List<String> fields) {
        StringBuilder line = new StringBuilder();
        for (int i = 0; i < fields.size(); i++) {
            if (i > 0) {
                line.append('\t');
            }
            line.append(fields.get(i));
        }
        return line.toString();
    }

    private static void writeRunInfo(
            File file,
            String inputFormat,
            int clusterDistance,
            String clusterRepresentative,
            int minGeneReads,
            int minSiteReads,
            double minPau,
            double minDeltaPau,
            double fdr,
            double wpiThreshold,
            double internalExonThreshold,
            double aleThreshold,
            double ipaThreshold,
            double tutrThreshold,
            String control,
            String caseCondition
    ) throws IOException {
        try (PrintWriter writer = new PrintWriter(new FileWriter(file))) {
            writer.println("APAURA long-read APA analysis");
            writer.println("input_format=" + inputFormat);
            writer.println("control=" + control);
            writer.println("case=" + caseCondition);
            writer.println("cluster_distance=" + clusterDistance);
            writer.println("cluster_representative=" + clusterRepresentative);
            writer.println("min_gene_reads_per_condition=" + minGeneReads);
            writer.println("min_site_reads_pooled=" + minSiteReads);
            writer.println("min_pau_pooled=" + minPau);
            writer.println("min_delta_pau=" + minDeltaPau);
            writer.println("fdr=" + fdr);
            writer.println("wpi_threshold=" + wpiThreshold);
            writer.println("min_internal_exon_delta=" + internalExonThreshold);
            writer.println("min_ale_switch=" + aleThreshold);
            writer.println("min_ipa_delta=" + ipaThreshold);
            writer.println("min_tutr_shift=" + tutrThreshold);
            writer.println(
                    "gene_level_test=tie-adjusted discrete two-sample ECDF maximum-distance test"
            );
        }
    }

    private static int count(SiteCluster site, String sampleName) {
        Integer value = site.counts.get(sampleName);
        return value == null ? 0 : value;
    }
}
