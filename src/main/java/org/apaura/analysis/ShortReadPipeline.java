package org.apaura.analysis;

import org.apaura.annotation.TranscriptAnnotation;
import org.apaura.io.SampleSheet;
import org.apaura.model.CoverageSegment;
import org.apaura.model.Sample;
import org.apaura.model.ShortReadGroup;
import org.apaura.model.ShortReadQuantification;
import org.apaura.model.ShortReadSite;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/** End-to-end conventional RNA-seq APA pipeline. */
public final class ShortReadPipeline {
    private ShortReadPipeline() {
    }

    public static void run(
            String sampleSheet,
            String gtf,
            String outDir,
            String control,
            String caseCondition,
            int minMapq,
            boolean keepDuplicates,
            String library,
            int anchorTolerance,
            int pasMergeDistance,
            int ipaWindow,
            boolean runIpa,
            int minSegmentLength,
            int minSampleEvidence,
            int minInformativeReplicates,
            boolean writePerSample,
            boolean writeDiagnostics,
            double minDeltaPau,
            double minIpaDelta,
            double fdr,
            double wpiThreshold
    ) throws Exception {
        File outputDirectory = ensureOutputDirectory(outDir);
        List<Sample> samples = SampleSheet.read(sampleSheet);
        validateConditions(samples, control, caseCondition);

        TranscriptAnnotation annotation = TranscriptAnnotation.load(gtf);
        writeArchitectureCatalog(
                new File(outputDirectory, "apa_architecture.tsv"),
                ArchitectureCatalog.build(
                        annotation.transcripts(),
                        anchorTolerance,
                        pasMergeDistance
                )
        );

        ShortReadGroupBuilder.BuildResult buildResult =
                ShortReadGroupBuilder.buildWithReport(
                        annotation.transcripts(),
                        anchorTolerance,
                        pasMergeDistance,
                        minSegmentLength
                );
        List<ShortReadGroup> groups = buildResult.groups;
        List<ShortReadIpa.Group> ipaGroups = runIpa
                ? ShortReadIpa.build(
                annotation.transcripts(),
                pasMergeDistance,
                ipaWindow,
                minSegmentLength
        )
                : Collections.<ShortReadIpa.Group>emptyList();

        if (groups.isEmpty() && ipaGroups.isEmpty()) {
            throw new IllegalArgumentException(
                    "No quantifiable terminal APA or IPA events were built from the GTF. "
                            + "Try lowering --min-segment-length, inspect --ipa-window, "
                            + "or inspect the annotation."
            );
        }

        System.err.println(
                "[APAURA] terminal-family reconciliation: initial_terminal_exon_groups="
                        + buildResult.initialAnchorGroups
                        + " reconciled_terminal_families="
                        + buildResult.reconciledFamilies
                        + " merged_terminal_exon_groups="
                        + buildResult.mergedAnchorGroups
                        + " skipped_nonquantifiable="
                        + buildResult.skippedNonQuantifiableFamilies
        );
        System.err.println("[APAURA] annotated tandem APA groups=" + groups.size());
        writeGroups(new File(outputDirectory, "pas.tsv"), groups);

        if (runIpa) {
            System.err.println(
                    "[APAURA] annotated IPA events="
                            + ipaGroups.size()
                            + " pas_merge_distance="
                            + pasMergeDistance
            );
            ShortReadIpa.writeAnnotation(
                    new File(outputDirectory, "ipa_annotation.tsv"),
                    ipaGroups
            );
        } else {
            System.err.println("[APAURA] IPA analysis=disabled");
            deleteIpaOutputs(outputDirectory);
        }

        QuantificationBundle quantifications = quantifySamples(
                samples,
                groups,
                ipaGroups,
                minMapq,
                keepDuplicates,
                library
        );

        if (writeDiagnostics) {
            writeSegments(
                    new File(outputDirectory, "pas_coverage.tsv"),
                    samples,
                    quantifications.apa
            );
        }

        if (writePerSample) {
            writeUsage(
                    new File(outputDirectory, "pau.tsv"),
                    samples,
                    quantifications.apa
            );
        }

        ShortReadDifferential.Results differential = ShortReadDifferential.run(
                samples,
                quantifications.apa,
                control,
                caseCondition,
                minSampleEvidence,
                minInformativeReplicates,
                minDeltaPau,
                fdr,
                wpiThreshold
        );
        writeSiteDiff(
                new File(outputDirectory, "apa_sites.tsv"),
                differential.sites
        );
        writeGroupSummary(
                new File(outputDirectory, "apa.tsv"),
                differential.groups
        );

        if (runIpa) {
            if (writePerSample) {
                ShortReadIpa.writeUsage(
                        new File(outputDirectory, "ipa_usage.tsv"),
                        samples,
                        quantifications.ipa
                );
            }

            ShortReadIpa.Results ipaDifferential = ShortReadIpa.differential(
                    samples,
                    quantifications.ipa,
                    control,
                    caseCondition,
                    minSampleEvidence,
                    minInformativeReplicates,
                    minIpaDelta,
                    fdr
            );
            ShortReadIpa.writeSites(
                    new File(outputDirectory, "ipa_sites.tsv"),
                    ipaDifferential.sites
            );
            ShortReadIpa.writeGenes(
                    new File(outputDirectory, "ipa.tsv"),
                    ipaDifferential.genes
            );
        }

        System.err.println(
                "[APAURA] short-read analysis complete: "
                        + outputDirectory.getAbsolutePath()
        );
    }

    private static void deleteIpaOutputs(File outputDirectory) throws IOException {
        String[] names = {
                "ipa_annotation.tsv",
                "ipa_usage.tsv",
                "ipa_sites.tsv",
                "ipa.tsv"
        };
        for (String name : names) {
            Files.deleteIfExists(new File(outputDirectory, name).toPath());
        }
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

    private static final class QuantificationBundle {
        final Map<String, List<ShortReadQuantification>> apa =
                new LinkedHashMap<String, List<ShortReadQuantification>>();
        final Map<String, List<ShortReadIpa.Quantification>> ipa =
                new LinkedHashMap<String, List<ShortReadIpa.Quantification>>();
    }

    private static QuantificationBundle quantifySamples(
            List<Sample> samples,
            List<ShortReadGroup> groups,
            List<ShortReadIpa.Group> ipaGroups,
            int minMapq,
            boolean keepDuplicates,
            String library
    ) throws Exception {
        QuantificationBundle quantifications = new QuantificationBundle();

        for (Sample sample : samples) {
            System.err.println(
                    "[APAURA] quantifying sample="
                            + sample.name
                            + " condition="
                            + sample.condition
            );
            ShortReadQuantifier.CombinedResult sampleResult =
                    ShortReadQuantifier.quantifyWithIpa(
                            sample.path,
                            groups,
                            ipaGroups,
                            minMapq,
                            keepDuplicates,
                            library
                    );
            quantifications.apa.put(sample.name, sampleResult.apa);
            quantifications.ipa.put(sample.name, sampleResult.ipa);
        }
        return quantifications;
    }

    private static void writeArchitectureCatalog(
            File file,
            List<ArchitectureCatalog.Entry> entries
    ) throws IOException {
        try (PrintWriter writer = writer(file)) {
            writer.println(
                    "gene_id\tgene_name\tbiotype\tchrom\tstrand\tannotation_architecture\t"
                            + "components\tn_transcripts\tn_terminal_branches\tn_annotated_pas\t"
                            + "n_tutr_branches\tn_ipa_compatible_ends"
            );

            for (ArchitectureCatalog.Entry entry : entries) {
                writer.println(
                        entry.geneId + "\t"
                                + entry.geneName + "\t"
                                + entry.biotype + "\t"
                                + entry.chrom + "\t"
                                + entry.strand.symbol() + "\t"
                                + entry.architectureClass + "\t"
                                + entry.components + "\t"
                                + entry.transcripts + "\t"
                                + entry.terminalBranches + "\t"
                                + entry.annotatedPas + "\t"
                                + entry.tutrBranches + "\t"
                                + entry.ipaCompatibleEnds
                );
            }
        }
    }

    private static void writeGroups(File file, List<ShortReadGroup> groups) throws IOException {
        try (PrintWriter writer = writer(file)) {
            writer.println(
                    "group_id\tgene_id\tgene_name\tbiotype\tchrom\tstrand\tscope\tanchor0\t"
                            + "n_sites\tpas_positions_1based\tsegment_lengths"
            );

            for (ShortReadGroup group : groups) {
                writer.println(
                        group.id + "\t"
                                + group.geneId + "\t"
                                + group.geneName + "\t"
                                + group.biotype + "\t"
                                + group.chrom + "\t"
                                + group.strand.symbol() + "\t"
                                + group.scope + "\t"
                                + group.anchor0 + "\t"
                                + group.sites.size() + "\t"
                                + commaSeparatedPositions(group) + "\t"
                                + commaSeparatedSegmentLengths(group)
                );
            }
        }
    }

    private static String commaSeparatedPositions(ShortReadGroup group) {
        StringBuilder values = new StringBuilder();
        for (int i = 0; i < group.sites.size(); i++) {
            if (i > 0) {
                values.append(',');
            }
            values.append(group.sites.get(i).pos0 + 1);
        }
        return values.toString();
    }

    private static String commaSeparatedSegmentLengths(ShortReadGroup group) {
        StringBuilder values = new StringBuilder();
        for (int i = 0; i < group.segments.size(); i++) {
            if (i > 0) {
                values.append(',');
            }
            values.append(group.segments.get(i).length());
        }
        return values.toString();
    }

    private static void writeSegments(
            File file,
            List<Sample> samples,
            Map<String, List<ShortReadQuantification>> quantifications
    ) throws IOException {
        try (PrintWriter writer = writer(file)) {
            writer.println(
                    "sample\tcondition\tgroup_id\tgene_id\tsegment_index\tsegment_label\tchrom\t"
                            + "start0\tend0\tlength\taligned_bases\traw_mean_depth\t"
                            + "fitted_mean_depth\tcommon_evidence_reads"
            );

            for (Sample sample : samples) {
                for (ShortReadQuantification quantification : quantifications.get(sample.name)) {
                    for (int segmentIndex = 0;
                         segmentIndex < quantification.group.segments.size();
                         segmentIndex++) {
                        CoverageSegment segment =
                                quantification.group.segments.get(segmentIndex);
                        writer.printf(
                                Locale.ROOT,
                                "%s\t%s\t%s\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%.8f\t%.8f\t%d%n",
                                sample.name,
                                sample.condition,
                                quantification.group.id,
                                quantification.group.geneId,
                                segmentIndex,
                                segment.label,
                                quantification.group.chrom,
                                segment.start0,
                                segment.end0,
                                segment.length(),
                                quantification.alignedBases[segmentIndex],
                                quantification.rawCoverage[segmentIndex],
                                quantification.fittedCoverage[segmentIndex],
                                quantification.commonEvidenceReads
                        );
                    }
                }
            }
        }
    }

    private static void writeUsage(
            File file,
            List<Sample> samples,
            Map<String, List<ShortReadQuantification>> quantifications
    ) throws IOException {
        try (PrintWriter writer = writer(file)) {
            writer.println(
                    "sample\tcondition\tgroup_id\tgene_id\tgene_name\tsite_id\tchrom\t"
                            + "pas_start0\tpas_end0\tpas_position_1based\tstrand\tsite_rank\t"
                            + "normalized_position\testimated_abundance\tpau\tpseudo_count\t"
                            + "total_evidence_reads"
            );

            for (Sample sample : samples) {
                for (ShortReadQuantification quantification : quantifications.get(sample.name)) {
                    for (int siteIndex = 0;
                         siteIndex < quantification.group.sites.size();
                         siteIndex++) {
                        ShortReadSite site = quantification.group.sites.get(siteIndex);
                        writer.printf(
                                Locale.ROOT,
                                "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%c\t%d\t%.6f\t%.8f\t%.8f\t%d\t%d%n",
                                sample.name,
                                sample.condition,
                                quantification.group.id,
                                quantification.group.geneId,
                                quantification.group.geneName,
                                site.id,
                                quantification.group.chrom,
                                site.pos0,
                                site.pos0 + 1,
                                site.pos0 + 1,
                                quantification.group.strand.symbol(),
                                site.rank,
                                site.normalizedPosition,
                                quantification.siteAbundance[siteIndex],
                                quantification.sitePau[siteIndex],
                                quantification.pseudoCounts[siteIndex],
                                quantification.commonEvidenceReads
                        );
                    }
                }
            }
        }
    }

    private static void writeSiteDiff(
            File file,
            List<ShortReadDifferential.SiteResult> results
    ) throws IOException {
        try (PrintWriter writer = writer(file)) {
            writer.println(
                    "group_id	gene_id	gene_name	biotype	scope	site_id	chrom	"
                            + "pas_position_1based	strand	site_rank	normalized_position	"
                            + "control_pau	case_pau	delta_pau	control_pau_sd	case_pau_sd	"
                            + "control_n_total	case_n_total	control_n_any_evidence	case_n_any_evidence	"
                            + "control_n_informative	case_n_informative	"
                            + "control_evidence	case_evidence	pvalue	fdr	method	status"
            );
            for (ShortReadDifferential.SiteResult result : results) {
                writer.println(
                        result.group.id + "	"
                                + result.group.geneId + "	"
                                + result.group.geneName + "	"
                                + result.group.biotype + "	"
                                + result.group.scope + "	"
                                + result.site.id + "	"
                                + result.group.chrom + "	"
                                + (result.site.pos0 + 1) + "	"
                                + result.group.strand.symbol() + "	"
                                + result.site.rank + "	"
                                + fmt(result.site.normalizedPosition, "%.6f") + "	"
                                + fmt(result.controlPau, "%.8f") + "	"
                                + fmt(result.casePau, "%.8f") + "	"
                                + fmt(result.delta, "%.8f") + "	"
                                + fmt(result.controlPauSd, "%.8f") + "	"
                                + fmt(result.casePauSd, "%.8f") + "	"
                                + result.controlTotalReplicates + "	"
                                + result.caseTotalReplicates + "	"
                                + result.controlAnyEvidence + "	"
                                + result.caseAnyEvidence + "	"
                                + result.controlInformativeReplicates + "	"
                                + result.caseInformativeReplicates + "	"
                                + result.controlEvidence + "	"
                                + result.caseEvidence + "	"
                                + fmt(result.p, "%.10g") + "	"
                                + fmt(result.q, "%.10g") + "	"
                                + result.method + "	"
                                + result.status
                );
            }
        }
    }

    private static void writeGroupSummary(
            File file,
            List<ShortReadDifferential.GroupResult> results
    ) throws IOException {
        try (PrintWriter writer = writer(file)) {
            writer.println(
                    "group_id	gene_id	gene_name	biotype	chrom	strand	scope	"
                            + "apa_architecture	n_sites	control_wpi	case_wpi	delta_wpi	"
                            + "control_wpi_sd	case_wpi_sd	apd	"
                            + "control_n_total	case_n_total	control_n_any_evidence	case_n_any_evidence	"
                            + "control_n_informative	case_n_informative	min_fdr	"
                            + "apa_change	observed_direction	apa_direction	apa_event	status"
            );
            for (ShortReadDifferential.GroupResult result : results) {
                writer.println(
                        result.group.id + "	"
                                + result.group.geneId + "	"
                                + result.group.geneName + "	"
                                + result.group.biotype + "	"
                                + result.group.chrom + "	"
                                + result.group.strand.symbol() + "	"
                                + result.group.scope + "	"
                                + result.apaArchitecture + "	"
                                + result.group.sites.size() + "	"
                                + fmt(result.controlWpi, "%.8f") + "	"
                                + fmt(result.caseWpi, "%.8f") + "	"
                                + fmt(result.deltaWpi, "%.8f") + "	"
                                + fmt(result.controlWpiSd, "%.8f") + "	"
                                + fmt(result.caseWpiSd, "%.8f") + "	"
                                + fmt(result.apd, "%.8f") + "	"
                                + result.controlTotalReplicates + "	"
                                + result.caseTotalReplicates + "	"
                                + result.controlAnyEvidence + "	"
                                + result.caseAnyEvidence + "	"
                                + result.controlInformativeReplicates + "	"
                                + result.caseInformativeReplicates + "	"
                                + fmt(result.minFdr, "%.10g") + "	"
                                + (result.apaChange ? "yes" : "no") + "	"
                                + result.observedDirection + "	"
                                + result.apaDirection + "	"
                                + result.apaEvent + "	"
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

    private static PrintWriter writer(File file) throws IOException {
        return new PrintWriter(new BufferedWriter(new FileWriter(file)));
    }
}