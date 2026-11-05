package org.apaura.analysis;

import org.apaura.annotation.TranscriptAnnotation;
import org.apaura.io.SampleSheet;
import org.apaura.output.ShortReadCanonicalAdapter;
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
            double minAleSwitch,
            double minIpaDelta,
            double fdr,
            double wpiThreshold,
            boolean discoverNovelPas,
            int novelPasWindow,
            double novelPasMinDepth,
            double novelPasMinDrop,
            int novelPasMinReplicates,
            int novelPasClusterDistance
    ) throws Exception {
        File outputDirectory = ensureOutputDirectory(outDir);
        File detailsDirectory = new File(outputDirectory, "details");
        if (writeDiagnostics || writePerSample || discoverNovelPas) {
            detailsDirectory = ensureSubdirectory(outputDirectory, "details");
        }
        List<Sample> samples = SampleSheet.read(sampleSheet);
        validateConditions(samples, control, caseCondition);

        TranscriptAnnotation annotation =
                TranscriptAnnotation.load(gtf);

        Map<String, List<Integer>> novelPasByFamily =
                Collections.<String, List<Integer>>emptyMap();

        int acceptedNovelPas = 0;

        if (discoverNovelPas) {

            System.err.println(
                    "[APAURA] novel TUTR PAS discovery=enabled"
                            + " method=coverage_breakpoint"
                            + " condition_blind=true"
            );

            ShortReadNovelPasCollector.Result discovery =
                    ShortReadNovelPasCollector.discover(
                            samples,
                            annotation.transcripts(),
                            anchorTolerance,
                            pasMergeDistance,
                            minMapq,
                            keepDuplicates,
                            library,
                            novelPasWindow,
                            novelPasMinDepth,
                            novelPasMinDrop,
                            novelPasMinReplicates,
                            novelPasClusterDistance,
                            minSegmentLength
                    );

            novelPasByFamily =
                    discovery.acceptedByFamily;

            for (List<Integer> positions
                    : novelPasByFamily.values()) {
                acceptedNovelPas += positions.size();
            }

            ShortReadNovelPasCollector.writeDiagnostics(
                    new File(
                            detailsDirectory,
                            "novel_pas_candidates.tsv"
                    ),
                    discovery
            );

            System.err.println(
                    "[APAURA] coverage-supported novel PAS accepted="
                            + acceptedNovelPas
                            + " candidate_rows="
                            + discovery.sites.size()
            );

        } else {

            Files.deleteIfExists(
                    new File(
                            detailsDirectory,
                            "novel_pas_candidates.tsv"
                    ).toPath()
            );

            System.err.println(
                    "[APAURA] novel TUTR PAS discovery=disabled"
            );
        }

        List<ArchitectureCatalog.Entry> architectureEntries =
                ArchitectureCatalog.build(
                        annotation.transcripts(),
                        anchorTolerance,
                        pasMergeDistance,
                        novelPasByFamily
                );

        if (writeDiagnostics) {
            writeArchitectureCatalog(
                    new File(
                            detailsDirectory,
                            "apa_architecture.tsv"
                    ),
                    architectureEntries
            );
        }

        ShortReadGroupBuilder.BuildResult buildResult =
                ShortReadGroupBuilder.buildWithReport(
                        annotation.transcripts(),
                        anchorTolerance,
                        pasMergeDistance,
                        minSegmentLength,
                        novelPasByFamily
                );
        List<ShortReadGroup> groups = buildResult.groups;
        List<ShortReadAle.Group> aleGroups = ShortReadAle.build(
                annotation.transcripts(),
                anchorTolerance,
                pasMergeDistance,
                minSegmentLength
        );
        List<ShortReadIpa.Group> ipaGroups = runIpa
                ? ShortReadIpa.build(
                annotation.transcripts(),
                pasMergeDistance,
                ipaWindow,
                minSegmentLength
        )
                : Collections.<ShortReadIpa.Group>emptyList();

        if (groups.isEmpty() && aleGroups.isEmpty() && ipaGroups.isEmpty()) {
            throw new IllegalArgumentException(
                    "No quantifiable terminal APA or IPA events were built from the annotation/discovery catalog. "
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
        System.err.println(
                "[APAURA] tandem APA groups="
                        + groups.size()
                        + " accepted_novel_pas="
                        + acceptedNovelPas
        );
        System.err.println("[APAURA] annotated ALE genes=" + aleGroups.size());
        if (writeDiagnostics) {
            writeGroups(new File(detailsDirectory, "pas.tsv"), groups);
            ShortReadAle.writeAnnotation(
                    new File(detailsDirectory, "ale_annotation.tsv"),
                    aleGroups
            );
        }

        if (runIpa) {
            System.err.println(
                    "[APAURA] annotated IPA events="
                            + ipaGroups.size()
                            + " pas_merge_distance="
                            + pasMergeDistance
            );
            if (writeDiagnostics) {
                ShortReadIpa.writeAnnotation(
                        new File(detailsDirectory, "ipa_annotation.tsv"),
                        ipaGroups
                );
            }
        } else {
            System.err.println("[APAURA] IPA analysis=disabled");
            deleteIpaOutputs(detailsDirectory);
        }

        QuantificationBundle quantifications = quantifySamples(
                samples,
                groups,
                ipaGroups,
                aleGroups,
                minMapq,
                keepDuplicates,
                library
        );

        if (writeDiagnostics) {
            writeSegments(
                    new File(detailsDirectory, "pas_coverage.tsv"),
                    samples,
                    quantifications.apa
            );
        }

        if (writePerSample) {
            writeUsage(
                    new File(detailsDirectory, "pau.tsv"),
                    samples,
                    quantifications.apa
            );
            ShortReadAle.writeUsage(
                    new File(detailsDirectory, "ale_usage.tsv"),
                    samples,
                    quantifications.ale
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
        writeSiteDiffCompact(
                new File(outputDirectory, "apa_sites.tsv"),
                differential.sites
        );
        if (writeDiagnostics) {
            writeSiteDiffDetailed(
                    new File(detailsDirectory, "apa_sites_full.tsv"),
                    differential.sites
            );
            writeGroupSummary(
                    new File(detailsDirectory, "apa_groups.tsv"),
                    differential.groups
            );
        }

        ShortReadAle.Results aleDifferential = ShortReadAle.differential(
                samples,
                quantifications.ale,
                control,
                caseCondition,
                minSampleEvidence,
                minInformativeReplicates,
                minDeltaPau,
                minAleSwitch,
                fdr,
                wpiThreshold
        );
        if (writeDiagnostics) {
            ShortReadAle.writeBranches(
                    new File(detailsDirectory, "ale_branches.tsv"),
                    aleDifferential.branches
            );
            ShortReadAle.writeGenes(
                    new File(detailsDirectory, "ale.tsv"),
                    aleDifferential.genes
            );
        }

        ShortReadIpa.Results ipaDifferential = new ShortReadIpa.Results();
        if (runIpa) {
            if (writePerSample) {
                ShortReadIpa.writeUsage(
                        new File(detailsDirectory, "ipa_usage.tsv"),
                        samples,
                        quantifications.ipa
                );
            }

            ipaDifferential = ShortReadIpa.differential(
                    samples,
                    quantifications.ipa,
                    control,
                    caseCondition,
                    minSampleEvidence,
                    minInformativeReplicates,
                    minIpaDelta,
                    fdr
            );
            if (writeDiagnostics) {
                ShortReadIpa.writeSites(
                        new File(detailsDirectory, "ipa_sites.tsv"),
                        ipaDifferential.sites
                );
                ShortReadIpa.writeGenes(
                        new File(detailsDirectory, "ipa.tsv"),
                        ipaDifferential.genes
                );
            }
        }

        writeEventsCompact(
                new File(outputDirectory, "apa_events.tsv"),
                differential.groups,
                aleDifferential.genes,
                ipaDifferential.sites
        );
        if (writeDiagnostics) {
            writeEventsDetailed(
                    new File(detailsDirectory, "apa_events_full.tsv"),
                    differential.groups,
                    aleDifferential.genes,
                    ipaDifferential.sites
            );
        }

        List<ShortReadIntegratedApa.Result> integrated =
                ShortReadIntegratedApa.integrate(
                        architectureEntries,
                        differential.groups,
                        aleDifferential.genes,
                        ipaDifferential.genes
                );
        ShortReadCanonicalAdapter.writeApa(
                new File(outputDirectory, "apa.tsv"),
                integrated,
                differential.groups,
                aleDifferential.genes,
                ipaDifferential.sites
        );
        if (writeDiagnostics) {
            writeIntegratedSummaryDetailed(
                    new File(detailsDirectory, "apa_full.tsv"),
                    integrated
            );
        }

        cleanupLegacyRootDetailOutputs(outputDirectory);
        cleanupUnrequestedDetails(detailsDirectory, writeDiagnostics, writePerSample);
        System.err.println("[APAURA] primary outputs: apa.tsv, apa_events.tsv, apa_sites.tsv");
        if (writeDiagnostics || writePerSample || discoverNovelPas) {
            System.err.println(
                    "[APAURA] optional details: "
                            + detailsDirectory.getAbsolutePath()
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

    private static File ensureSubdirectory(File parent, String name) throws IOException {
        File directory = new File(parent, name);
        if (!directory.exists() && !directory.mkdirs()) {
            throw new IOException("Cannot create output directory " + directory);
        }
        return directory;
    }

    private static void cleanupLegacyRootDetailOutputs(File outputDirectory) throws IOException {
        String[] names = {
                "apa_architecture.tsv",
                "pas.tsv",
                "apa_groups.tsv",
                "ale_annotation.tsv",
                "ale_branches.tsv",
                "ale.tsv",
                "ipa_annotation.tsv",
                "ipa_sites.tsv",
                "ipa.tsv",
                "pau.tsv",
                "ale_usage.tsv",
                "ipa_usage.tsv",
                "pas_coverage.tsv"
        };
        for (String name : names) {
            Files.deleteIfExists(new File(outputDirectory, name).toPath());
        }
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
        final Map<String, List<ShortReadAle.Quantification>> ale =
                new LinkedHashMap<String, List<ShortReadAle.Quantification>>();
    }

    private static QuantificationBundle quantifySamples(
            List<Sample> samples,
            List<ShortReadGroup> groups,
            List<ShortReadIpa.Group> ipaGroups,
            List<ShortReadAle.Group> aleGroups,
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
                    ShortReadQuantifier.quantifyWithIpaAndAle(
                            sample.path,
                            groups,
                            ipaGroups,
                            aleGroups,
                            minMapq,
                            keepDuplicates,
                            library
                    );
            quantifications.apa.put(sample.name, sampleResult.apa);
            quantifications.ipa.put(sample.name, sampleResult.ipa);
            quantifications.ale.put(sample.name, sampleResult.ale);
        }
        return quantifications;
    }

    private static void writeArchitectureCatalog(
            File file,
            List<ArchitectureCatalog.Entry> entries
    ) throws IOException {
        try (PrintWriter writer = writer(file)) {
            writer.println(
                    "gene_id\tgene_name\tbiotype\tchrom\tgene_start_1based\tgene_end_1based\t"
                            + "strand\tannotation_architecture\tcomponents\tn_transcripts\t"
                            + "n_terminal_branches\tn_annotated_pas\tn_terminal_pas\t"
                            + "n_tutr_branches\tn_ipa_compatible_ends"
            );

            for (ArchitectureCatalog.Entry entry : entries) {
                writer.println(
                        entry.geneId + "\t"
                                + entry.geneName + "\t"
                                + entry.biotype + "\t"
                                + entry.chrom + "\t"
                                + entry.geneStart1based + "\t"
                                + entry.geneEnd1based + "\t"
                                + entry.strand.symbol() + "\t"
                                + entry.architectureClass + "\t"
                                + entry.components + "\t"
                                + entry.transcripts + "\t"
                                + entry.terminalBranches + "\t"
                                + entry.annotatedPas + "\t"
                                + entry.terminalPas + "\t"
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

    private static void writeSiteDiffCompact(
            File file,
            List<ShortReadDifferential.SiteResult> results
    ) throws IOException {
        try (PrintWriter writer = writer(file)) {
            writer.println(
                    "group_id\tgene_id\tgene_name\tchrom\tstrand\tsite_id\t"
                            + "pas_position_1based\tsite_rank\tcontrol_pau\tcase_pau\tdelta_pau\t"
                            + "pvalue\tfdr\tstatus"
            );
            for (ShortReadDifferential.SiteResult result : results) {
                writer.println(
                        result.group.id + "\t"
                                + result.group.geneId + "\t"
                                + result.group.geneName + "\t"
                                + result.group.chrom + "\t"
                                + result.group.strand.symbol() + "\t"
                                + result.site.id + "\t"
                                + (result.site.pos0 + 1) + "\t"
                                + result.site.rank + "\t"
                                + fmt(result.controlPau, "%.8f") + "\t"
                                + fmt(result.casePau, "%.8f") + "\t"
                                + fmt(result.delta, "%.8f") + "\t"
                                + fmt(result.p, "%.10g") + "\t"
                                + fmt(result.q, "%.10g") + "\t"
                                + result.status
                );
            }
        }
    }

    private static void writeSiteDiffDetailed(
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
                    "group_id\tgene_id\tgene_name\tbiotype\tchrom\tstrand\tscope\t"
                            + "apa_architecture\tn_sites\tcontrol_wpi\tcase_wpi\tdelta_wpi\t"
                            + "control_wpi_sd\tcase_wpi_sd\tapd\tjsd\t"
                            + "control_n_total\tcase_n_total\tcontrol_n_any_evidence\tcase_n_any_evidence\t"
                            + "control_n_informative\tcase_n_informative\tmin_pvalue\tmin_fdr\t"
                            + "apa_change\tobserved_direction\tapa_direction\tapa_event\tstatus"
            );
            for (ShortReadDifferential.GroupResult result : results) {
                writer.println(
                        result.group.id + "\t"
                                + result.group.geneId + "\t"
                                + result.group.geneName + "\t"
                                + result.group.biotype + "\t"
                                + result.group.chrom + "\t"
                                + result.group.strand.symbol() + "\t"
                                + result.group.scope + "\t"
                                + result.apaArchitecture + "\t"
                                + result.group.sites.size() + "\t"
                                + fmt(result.controlWpi, "%.8f") + "\t"
                                + fmt(result.caseWpi, "%.8f") + "\t"
                                + fmt(result.deltaWpi, "%.8f") + "\t"
                                + fmt(result.controlWpiSd, "%.8f") + "\t"
                                + fmt(result.caseWpiSd, "%.8f") + "\t"
                                + fmt(result.apd, "%.8f") + "\t"
                                + fmt(result.jsd, "%.8f") + "\t"
                                + result.controlTotalReplicates + "\t"
                                + result.caseTotalReplicates + "\t"
                                + result.controlAnyEvidence + "\t"
                                + result.caseAnyEvidence + "\t"
                                + result.controlInformativeReplicates + "\t"
                                + result.caseInformativeReplicates + "\t"
                                + fmt(result.minPvalue, "%.10g") + "\t"
                                + fmt(result.minFdr, "%.10g") + "\t"
                                + (result.apaChange ? "yes" : "no") + "\t"
                                + result.observedDirection + "\t"
                                + result.apaDirection + "\t"
                                + result.apaEvent + "\t"
                                + result.status
                );
            }
        }
    }

    private static void writeEventsCompact(
            File file,
            List<ShortReadDifferential.GroupResult> tutrGroups,
            List<ShortReadAle.GeneResult> aleGenes,
            List<ShortReadIpa.SiteResult> ipaSites
    ) throws IOException {
        try (PrintWriter writer = writer(file)) {
            writer.println(
                    "event_id\tgene_id\tgene_name\tbiotype\tchrom\tstrand\tevent_type\t"
                            + "start_1based\tend_1based\tcontrol_usage\tcase_usage\tdelta_usage\t"
                            + "effect_size\tpvalue\tfdr\tdirection\tstatus"
            );

            for (ShortReadDifferential.GroupResult result : tutrGroups) {
                int start1 = Integer.MAX_VALUE;
                int end1 = Integer.MIN_VALUE;
                for (ShortReadSite site : result.group.sites) {
                    int position1 = site.pos0 + 1;
                    start1 = Math.min(start1, position1);
                    end1 = Math.max(end1, position1);
                }
                writer.println(
                        result.group.id + "\t" + result.group.geneId + "\t"
                                + result.group.geneName + "\t" + result.group.biotype + "\t"
                                + result.group.chrom + "\t" + result.group.strand.symbol() + "\tTUTR\t"
                                + coordinate(start1) + "\t" + coordinate(end1) + "\t"
                                + fmt(result.controlWpi, "%.8f") + "\t"
                                + fmt(result.caseWpi, "%.8f") + "\t"
                                + fmt(result.deltaWpi, "%.8f") + "\t"
                                + fmt(result.apd, "%.8f") + "\t"
                                + fmt(result.minPvalue, "%.10g") + "\t"
                                + fmt(result.minFdr, "%.10g") + "\t"
                                + userDirection(result.apaDirection) + "\t" + result.status
                );
            }

            for (ShortReadAle.GeneResult result : aleGenes) {
                List<Integer> positions = alePositions1(result.group);
                int start1 = positions.isEmpty() ? Integer.MAX_VALUE : positions.get(0);
                int end1 = positions.isEmpty() ? Integer.MIN_VALUE : positions.get(positions.size() - 1);
                writer.println(
                        result.group.id + "\t" + result.group.geneId + "\t"
                                + result.group.geneName + "\t" + result.group.biotype + "\t"
                                + result.group.chrom + "\t" + result.group.strand.symbol() + "\tALE\t"
                                + coordinate(start1) + "\t" + coordinate(end1) + "\t"
                                + fmt(result.controlWpi, "%.8f") + "\t"
                                + fmt(result.caseWpi, "%.8f") + "\t"
                                + fmt(result.deltaWpi, "%.8f") + "\t"
                                + fmt(result.switchIndex, "%.8f") + "\t"
                                + fmt(result.minPvalue, "%.10g") + "\t"
                                + fmt(result.minFdr, "%.10g") + "\t"
                                + userDirection(result.aleDirection) + "\t" + result.status
                );
            }

            for (ShortReadIpa.SiteResult result : ipaSites) {
                ShortReadIpa.Group group = result.group;
                writer.println(
                        group.id + "\t" + group.geneId + "\t" + group.geneName + "\t"
                                + group.biotype + "\t" + group.chrom + "\t" + group.strand.symbol() + "\tIPA\t"
                                + (group.intronStart0 + 1) + "\t" + group.intronEnd0 + "\t"
                                + fmt(result.controlUsage, "%.8f") + "\t"
                                + fmt(result.caseUsage, "%.8f") + "\t"
                                + fmt(result.delta, "%.8f") + "\t"
                                + fmt(Double.isFinite(result.delta) ? Math.abs(result.delta) : Double.NaN, "%.8f")
                                + "\t"
                                + fmt(result.p, "%.10g") + "\t" + fmt(result.q, "%.10g") + "\t"
                                + userDirection(result.direction) + "\t" + result.status
                );
            }
        }
    }

    private static void writeEventsDetailed(
            File file,
            List<ShortReadDifferential.GroupResult> tutrGroups,
            List<ShortReadAle.GeneResult> aleGenes,
            List<ShortReadIpa.SiteResult> ipaSites
    ) throws IOException {
        try (PrintWriter writer = writer(file)) {
            writer.println(
                    "event_id\tgene_id\tgene_name\tbiotype\tchrom\tstrand\tevent_type\tscope\t"
                            + "start_1based\tend_1based\trepresentative_position_1based\t"
                            + "n_sites\tpositions_1based\tcontrol_usage\tcase_usage\tdelta_usage\t"
                            + "effect_size\tpvalue\tfdr\tobserved_direction\tdirection\tmethod\tstatus"
            );

            for (ShortReadDifferential.GroupResult result : tutrGroups) {
                int start1 = Integer.MAX_VALUE;
                int end1 = Integer.MIN_VALUE;
                for (ShortReadSite site : result.group.sites) {
                    int position1 = site.pos0 + 1;
                    start1 = Math.min(start1, position1);
                    end1 = Math.max(end1, position1);
                }
                writer.println(
                        result.group.id + "\t"
                                + result.group.geneId + "\t"
                                + result.group.geneName + "\t"
                                + result.group.biotype + "\t"
                                + result.group.chrom + "\t"
                                + result.group.strand.symbol() + "\t"
                                + "TUTR\t"
                                + result.group.scope + "\t"
                                + coordinate(start1) + "\t"
                                + coordinate(end1) + "\t"
                                + "NA\t"
                                + result.group.sites.size() + "\t"
                                + commaSeparatedPositions(result.group) + "\t"
                                + fmt(result.controlWpi, "%.8f") + "\t"
                                + fmt(result.caseWpi, "%.8f") + "\t"
                                + fmt(result.deltaWpi, "%.8f") + "\t"
                                + fmt(result.apd, "%.8f") + "\t"
                                + fmt(result.minPvalue, "%.10g") + "\t"
                                + fmt(result.minFdr, "%.10g") + "\t"
                                + result.observedDirection + "\t"
                                + result.apaDirection + "\t"
                                + "min_of_PAS_tests\t"
                                + result.status
                );
            }

            for (ShortReadAle.GeneResult result : aleGenes) {
                List<Integer> positions = alePositions1(result.group);
                int start1 = positions.isEmpty() ? Integer.MAX_VALUE : positions.get(0);
                int end1 = positions.isEmpty() ? Integer.MIN_VALUE : positions.get(positions.size() - 1);
                writer.println(
                        result.group.id + "\t"
                                + result.group.geneId + "\t"
                                + result.group.geneName + "\t"
                                + result.group.biotype + "\t"
                                + result.group.chrom + "\t"
                                + result.group.strand.symbol() + "\t"
                                + "ALE\tALTERNATIVE_LAST_EXON\t"
                                + coordinate(start1) + "\t"
                                + coordinate(end1) + "\t"
                                + "NA\t"
                                + positions.size() + "\t"
                                + joinIntegers(positions) + "\t"
                                + fmt(result.controlWpi, "%.8f") + "\t"
                                + fmt(result.caseWpi, "%.8f") + "\t"
                                + fmt(result.deltaWpi, "%.8f") + "\t"
                                + fmt(result.switchIndex, "%.8f") + "\t"
                                + fmt(result.minPvalue, "%.10g") + "\t"
                                + fmt(result.minFdr, "%.10g") + "\t"
                                + result.observedDirection + "\t"
                                + result.aleDirection + "\t"
                                + "min_of_ALE_branch_tests\t"
                                + result.status
                );
            }

            for (ShortReadIpa.SiteResult result : ipaSites) {
                ShortReadIpa.Group group = result.group;
                List<Integer> positions = new java.util.ArrayList<Integer>();
                for (Integer position0 : group.annotatedPas0) {
                    positions.add(position0 + 1);
                }
                java.util.Collections.sort(positions);
                writer.println(
                        group.id + "\t"
                                + group.geneId + "\t"
                                + group.geneName + "\t"
                                + group.biotype + "\t"
                                + group.chrom + "\t"
                                + group.strand.symbol() + "\t"
                                + "IPA\tINTRONIC_POLYADENYLATION\t"
                                + (group.intronStart0 + 1) + "\t"
                                + group.intronEnd0 + "\t"
                                + (group.representativePas0 + 1) + "\t"
                                + positions.size() + "\t"
                                + joinIntegers(positions) + "\t"
                                + fmt(result.controlUsage, "%.8f") + "\t"
                                + fmt(result.caseUsage, "%.8f") + "\t"
                                + fmt(result.delta, "%.8f") + "\t"
                                + fmt(Double.isFinite(result.delta) ? Math.abs(result.delta) : Double.NaN, "%.8f") + "\t"
                                + fmt(result.p, "%.10g") + "\t"
                                + fmt(result.q, "%.10g") + "\t"
                                + observedIpaDirection(result.delta) + "\t"
                                + result.direction + "\t"
                                + result.method + "\t"
                                + result.status
                );
            }
        }
    }

    private static String observedIpaDirection(double delta) {
        if (!Double.isFinite(delta)) return "not_applicable";
        if (delta > 0.0) return "gain";
        if (delta < 0.0) return "loss";
        return "no_change";
    }

    private static List<Integer> alePositions1(ShortReadAle.Group group) {
        Set<Integer> positions = new LinkedHashSet<Integer>();
        for (ShortReadAle.Branch branch : group.branches) {
            for (Integer position0 : branch.annotatedPas0) {
                positions.add(position0 + 1);
            }
        }
        List<Integer> ordered = new java.util.ArrayList<Integer>(positions);
        java.util.Collections.sort(ordered);
        return ordered;
    }

    private static String joinIntegers(List<Integer> values) {
        StringBuilder out = new StringBuilder();
        for (Integer value : values) {
            if (out.length() > 0) out.append(',');
            out.append(value);
        }
        return out.toString();
    }

    private static String coordinate(int value) {
        return value == Integer.MAX_VALUE || value == Integer.MIN_VALUE
                ? "NA"
                : Integer.toString(value);
    }


    private static void writeIntegratedSummaryDetailed(
            File file,
            List<ShortReadIntegratedApa.Result> results
    ) throws IOException {
        try (PrintWriter writer = writer(file)) {
            writer.println(
                    "gene_id\tgene_name\tbiotype\tchrom\tgene_start_1based\tgene_end_1based\tstrand\t"
                            + "n_pas\tn_terminal_pas\tn_transcripts\tn_terminal_branches\tn_tutr_branches\t"
                            + "n_tutr_groups\tn_tutr_changed_groups\tn_ipa_events\tn_ipa_changed_events\t"
                            + "apa_architecture\tarchitecture_components\t"
                            + "tutr_representative_event\ttutr_change\ttutr_observed_direction\ttutr_direction\t"
                            + "delta_tutr_wpi\ttutr_shift_magnitude\ttutr_jsd\ttutr_min_pvalue\ttutr_min_fdr\ttutr_status\t"
                            + "ale_change\tale_observed_direction\tale_direction\tdelta_ale_wpi\tale_switch_index\t"
                            + "ale_min_pvalue\tale_min_fdr\tale_status\t"
                            + "ipa_representative_event\tipa_change\tipa_direction\tdelta_ipa_usage\t"
                            + "ipa_max_abs_delta\tipa_min_pvalue\tipa_min_fdr\tipa_status\t"
                            + "changed_mechanisms\tmechanism_directions\tapa_change\t"
                            + "apa_direction\tapa_event\tmin_mechanism_pvalue\tmin_mechanism_fdr\tstatus"
            );
            for (ShortReadIntegratedApa.Result result : results) {
                ArchitectureCatalog.Entry a = result.architecture;
                writer.println(
                        a.geneId + "\t"
                                + a.geneName + "\t"
                                + a.biotype + "\t"
                                + a.chrom + "\t"
                                + a.geneStart1based + "\t"
                                + a.geneEnd1based + "\t"
                                + a.strand.symbol() + "\t"
                                + a.annotatedPas + "\t"
                                + a.terminalPas + "\t"
                                + a.transcripts + "\t"
                                + a.terminalBranches + "\t"
                                + a.tutrBranches + "\t"
                                + result.nTutrGroups + "\t"
                                + result.nTutrChangedGroups + "\t"
                                + result.nIpaEvents + "\t"
                                + result.nIpaChangedEvents + "\t"
                                + a.architectureClass + "\t"
                                + a.components + "\t"
                                + result.tutrRepresentativeGroup + "\t"
                                + (result.tutrChange ? "yes" : "no") + "\t"
                                + result.tutrObservedDirection + "\t"
                                + result.tutrDirection + "\t"
                                + fmt(result.tutrDeltaWpi, "%.8f") + "\t"
                                + fmt(result.tutrShiftMagnitude, "%.8f") + "\t"
                                + fmt(result.tutrJsd, "%.8f") + "\t"
                                + fmt(result.tutrMinPvalue, "%.10g") + "\t"
                                + fmt(result.tutrMinFdr, "%.10g") + "\t"
                                + result.tutrStatus + "\t"
                                + (result.aleChange ? "yes" : "no") + "\t"
                                + result.aleObservedDirection + "\t"
                                + result.aleDirection + "\t"
                                + fmt(result.aleDeltaWpi, "%.8f") + "\t"
                                + fmt(result.aleSwitchIndex, "%.8f") + "\t"
                                + fmt(result.aleMinPvalue, "%.10g") + "\t"
                                + fmt(result.aleMinFdr, "%.10g") + "\t"
                                + result.aleStatus + "\t"
                                + result.ipaRepresentativeEvent + "\t"
                                + (result.ipaChange ? "yes" : "no") + "\t"
                                + result.ipaDirection + "\t"
                                + fmt(result.ipaDeltaUsage, "%.8f") + "\t"
                                + fmt(result.ipaMaxAbsDelta, "%.8f") + "\t"
                                + fmt(result.ipaMinPvalue, "%.10g") + "\t"
                                + fmt(result.ipaMinFdr, "%.10g") + "\t"
                                + result.ipaStatus + "\t"
                                + result.changedMechanisms + "\t"
                                + result.mechanismDirections + "\t"
                                + (result.apaChange ? "yes" : "no") + "\t"
                                + result.apaDirection + "\t"
                                + result.apaEvent + "\t"
                                + fmt(result.minPvalue, "%.10g") + "\t"
                                + fmt(result.minFdr, "%.10g") + "\t"
                                + result.status
                );
            }
        }
    }

    private static String userDirection(String direction) {
        if (direction == null) return "not_applicable";
        if ("redistribution".equals(direction) || "switch".equals(direction)) {
            return "redistribution";
        }
        return direction;
    }

    private static void cleanupUnrequestedDetails(
            File detailsDirectory,
            boolean writeDiagnostics,
            boolean writePerSample
    ) throws IOException {
        if (!detailsDirectory.exists()) return;

        if (!writeDiagnostics) {
            String[] diagnostics = {
                    "apa_architecture.tsv", "pas.tsv", "apa_groups.tsv",
                    "ale_annotation.tsv", "ale_branches.tsv", "ale.tsv",
                    "ipa_annotation.tsv", "ipa_sites.tsv", "ipa.tsv",
                    "pas_coverage.tsv", "apa_full.tsv", "apa_events_full.tsv",
                    "apa_sites_full.tsv"
            };
            for (String name : diagnostics) {
                Files.deleteIfExists(new File(detailsDirectory, name).toPath());
            }
        }

        if (!writePerSample) {
            String[] perSample = {"pau.tsv", "ale_usage.tsv", "ipa_usage.tsv"};
            for (String name : perSample) {
                Files.deleteIfExists(new File(detailsDirectory, name).toPath());
            }
        }

        String[] remaining = detailsDirectory.list();
        if (remaining != null && remaining.length == 0) {
            Files.deleteIfExists(detailsDirectory.toPath());
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