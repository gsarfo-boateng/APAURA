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
            int minSegmentLength,
            int minSampleEvidence,
            double minDeltaPau,
            double fdr,
            double wpiThreshold
    ) throws Exception {
        File outputDirectory = ensureOutputDirectory(outDir);
        List<Sample> samples = SampleSheet.read(sampleSheet);
        validateConditions(samples, control, caseCondition);

        TranscriptAnnotation annotation = TranscriptAnnotation.load(gtf);
        writeArchitectureCatalog(
                new File(outputDirectory, "apa_architecture.tsv"),
                ArchitectureCatalog.build(annotation.transcripts())
        );

        List<ShortReadGroup> groups = ShortReadGroupBuilder.build(
                annotation.transcripts(),
                anchorTolerance,
                pasMergeDistance,
                minSegmentLength
        );
        if (groups.isEmpty()) {
            throw new IllegalArgumentException(
                    "No quantifiable tandem terminal-exon APA groups were built from the GTF. "
                            + "Try lowering --min-segment-length or --pas-merge-distance, "
                            + "or inspect the annotation."
            );
        }

        System.err.println("[APAURA] annotated tandem APA groups=" + groups.size());
        writeGroups(new File(outputDirectory, "pas.tsv"), groups);

        Map<String, List<ShortReadQuantification>> quantifications = quantifySamples(
                samples,
                groups,
                minMapq,
                keepDuplicates,
                library
        );

        writeSegments(
                new File(outputDirectory, "pas_coverage.tsv"),
                samples,
                quantifications
        );
        writeUsage(
                new File(outputDirectory, "pau.tsv"),
                samples,
                quantifications
        );

        ShortReadDifferential.Results differential = ShortReadDifferential.run(
                samples,
                quantifications,
                control,
                caseCondition,
                minSampleEvidence,
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

        System.err.println(
                "[APAURA] short-read analysis complete: "
                        + outputDirectory.getAbsolutePath()
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

    private static Map<String, List<ShortReadQuantification>> quantifySamples(
            List<Sample> samples,
            List<ShortReadGroup> groups,
            int minMapq,
            boolean keepDuplicates,
            String library
    ) throws Exception {
        Map<String, List<ShortReadQuantification>> quantifications =
                new LinkedHashMap<String, List<ShortReadQuantification>>();

        for (Sample sample : samples) {
            System.err.println(
                    "[APAURA] quantifying sample="
                            + sample.name
                            + " condition="
                            + sample.condition
            );
            quantifications.put(
                    sample.name,
                    ShortReadQuantifier.quantify(
                            sample.path,
                            groups,
                            minMapq,
                            keepDuplicates,
                            library
                    )
            );
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
                    "group_id\tgene_id\tgene_name\tbiotype\tscope\tsite_id\tchrom\t"
                            + "pas_start0\tpas_end0\tpas_position_1based\tstrand\tsite_rank\t"
                            + "normalized_position\tcontrol_pau\tcase_pau\tdelta_pau\t"
                            + "control_evidence\tcase_evidence\tpvalue\tfdr\tdispersion_rho\tmethod"
            );

            for (ShortReadDifferential.SiteResult result : results) {
                writer.printf(
                        Locale.ROOT,
                        "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%c\t%d\t"
                                + "%.6f\t%.8f\t%.8f\t%.8f\t%d\t%d\t%.10g\t%.10g\t%s\t%s%n",
                        result.group.id,
                        result.group.geneId,
                        result.group.geneName,
                        result.group.biotype,
                        result.group.scope,
                        result.site.id,
                        result.group.chrom,
                        result.site.pos0,
                        result.site.pos0 + 1,
                        result.site.pos0 + 1,
                        result.group.strand.symbol(),
                        result.site.rank,
                        result.site.normalizedPosition,
                        result.controlPau,
                        result.casePau,
                        result.delta,
                        result.controlEvidence,
                        result.caseEvidence,
                        result.p,
                        result.q,
                        Double.isNaN(result.rho)
                                ? "NA"
                                : String.format(Locale.ROOT, "%.8f", result.rho),
                        result.method
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
                    "group_id\tgene_id\tgene_name\tbiotype\tchrom\tstrand\tscope\t"
                            + "apa_architecture\tapa_change\tapa_direction\tapa_event\tn_sites\t"
                            + "control_wpi\tcase_wpi\tdelta_wpi\tapd\tmin_fdr\tstatus"
            );

            for (ShortReadDifferential.GroupResult result : results) {
                writer.printf(
                        Locale.ROOT,
                        "%s\t%s\t%s\t%s\t%s\t%c\t%s\t%s\t%s\t%s\t%s\t%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.10g\t%s%n",
                        result.group.id,
                        result.group.geneId,
                        result.group.geneName,
                        result.group.biotype,
                        result.group.chrom,
                        result.group.strand.symbol(),
                        result.group.scope,
                        result.apaArchitecture,
                        result.apaChange ? "yes" : "no",
                        result.apaDirection,
                        result.apaEvent,
                        result.group.sites.size(),
                        result.controlWpi,
                        result.caseWpi,
                        result.deltaWpi,
                        result.apd,
                        result.minFdr,
                        result.status
                );
            }
        }
    }

    private static PrintWriter writer(File file) throws IOException {
        return new PrintWriter(new BufferedWriter(new FileWriter(file)));
    }
}
