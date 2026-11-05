package org.apaura.sc;

import org.apaura.analysis.ApaUsageMetrics;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Adds APAURA-native usage summaries to preserved sc3UTR differential output.
 *
 * <p>The embedded sc3UTR engine remains responsible for molecule counting,
 * biological-sample replicate handling, formal tests, interaction analysis,
 * and its legacy outputs. This harmonizer independently reconstructs
 * biological-sample PAS-usage vectors from the same input tables and reports
 * mean PAU, total-variation distance (APD), Jensen-Shannon divergence, and
 * weighted-position index (WPI) without treating cells as replicates.</p>
 */
final class ScApaUsageHarmonizer {

    static final String GENE_OUTPUT = "apaura_usage.tsv";
    static final String SITE_OUTPUT = "apaura_site_usage.tsv";

    private ScApaUsageHarmonizer() {
        // Utility class.
    }

    static void run(String[] args) throws IOException {
        Path sitesPath = requiredPath(args, "--sites");
        Path samplesPath = requiredPath(args, "--samples");
        Path countsPath = requiredPath(args, "--counts");
        Path outDir = requiredPath(args, "--out");
        String control = requiredArg(args, "--control");
        String treated = requiredArg(args, "--treated");

        int minSampleCount = intArg(args, "--min-sample-count", 5);
        double legacyMinDelta = doubleArg(args, "--min-delta", 0.10);
        double minTvd = doubleArg(args, "--min-tvd", legacyMinDelta);
        double wpiThreshold = doubleArg(
                args,
                "--wpi-threshold",
                legacyMinDelta
        );
        double fdr = doubleArg(args, "--fdr", 0.05);

        if (minSampleCount < 0) {
            throw new IllegalArgumentException(
                    "--min-sample-count must be >= 0"
            );
        }
        requireProbabilityThreshold("--min-tvd", minTvd);
        requireProbabilityThreshold("--wpi-threshold", wpiThreshold);
        requireProbabilityThreshold("--fdr", fdr);

        Path legacyResults = outDir.resolve("all_results.tsv");
        if (!Files.isRegularFile(legacyResults)) {
            throw new IOException(
                    "Single-cell comparison output is missing: "
                            + legacyResults
            );
        }

        Map<String, Site> sites = loadSites(sitesPath);
        Map<String, List<Site>> genes = groupSites(sites.values());
        Map<String, Sample> samples = loadSamples(samplesPath);
        validateConditions(samples.values(), control, treated);

        Grouping grouping = new Grouping(
                groupByArg(args),
                stringArg(args, "--group-separator", "::"),
                hasFlag(args, "--label-group-columns")
        );

        Path cellsPath = optionalPath(args, "--cells");
        Map<CountKey, Map<String, Long>> counts = loadCounts(
                countsPath,
                sites,
                cellsPath,
                grouping
        );

        TsvTable legacy = TsvTable.read(legacyResults);
        legacy.require(
                "group",
                "gene_id",
                "gene_symbol",
                "chrom",
                "strand",
                "n_sites",
                "analysis_mode",
                "n_control_samples",
                "n_treated_samples",
                "control_total",
                "treated_total",
                "control_utr_index",
                "treated_utr_index",
                "delta_utr_index",
                "pvalue",
                "qvalue",
                "distribution_pvalue",
                "distribution_qvalue",
                "test",
                "distribution_test",
                "classification"
        );

        Files.createDirectories(outDir);
        Path geneOutput = outDir.resolve(GENE_OUTPUT);
        Path siteOutput = outDir.resolve(SITE_OUTPUT);

        try (BufferedWriter geneWriter = Files.newBufferedWriter(
                geneOutput,
                StandardCharsets.UTF_8
        ); BufferedWriter siteWriter = Files.newBufferedWriter(
                siteOutput,
                StandardCharsets.UTF_8
        )) {
            writeGeneHeader(geneWriter);
            writeSiteHeader(siteWriter);

            for (String[] row : legacy.rows) {
                String group = legacy.get(row, "group");
                String geneId = legacy.get(row, "gene_id");
                List<Site> geneSites = genes.get(geneId);

                if (geneSites == null || geneSites.size() < 2) {
                    continue;
                }

                UsageSummary summary = summarize(
                        group,
                        geneSites,
                        samples.values(),
                        counts,
                        control,
                        treated,
                        minSampleCount
                );
                validateAgainstLegacy(
                        legacy,
                        row,
                        geneSites,
                        summary
                );

                String analysisMode = legacy.get(row, "analysis_mode");
                String primaryTest = legacy.get(row, "test");
                String distributionTest = legacy.get(
                        row,
                        "distribution_test"
                );
                double primaryP = parseDouble(
                        legacy.get(row, "pvalue")
                );
                double primaryQ = parseDouble(
                        legacy.get(row, "qvalue")
                );
                double distributionP = parseDouble(
                        legacy.get(row, "distribution_pvalue")
                );
                double distributionQ = parseDouble(
                        legacy.get(row, "distribution_qvalue")
                );

                boolean replicated =
                        "REPLICATED_DIFFERENTIAL".equals(analysisMode);
                boolean directionalSupport = replicated
                        && finiteAtMost(primaryQ, fdr);
                boolean distributionSupport = replicated
                        && geneSites.size() > 2
                        && "SAMPLE_LABEL_PERMUTATION_JS".equals(
                                distributionTest
                        )
                        && finiteAtMost(distributionQ, fdr);
                boolean formalSupport =
                        directionalSupport || distributionSupport;
                boolean observedChange = Double.isFinite(summary.apd)
                        && summary.apd >= minTvd;
                boolean apaChange = formalSupport && observedChange;

                String observedDirection =
                        ApaUsageMetrics.observedDirection(
                                summary.deltaWpi,
                                wpiThreshold
                        );
                String apaDirection =
                        ApaUsageMetrics.supportedDirection(
                                observedChange,
                                summary.deltaWpi,
                                wpiThreshold
                        );

                String evidenceStatus = evidenceStatus(
                        analysisMode,
                        primaryTest,
                        formalSupport
                );
                String apaClassification = classification(
                        apaChange,
                        apaDirection,
                        false
                );
                String observedClassification = classification(
                        observedChange,
                        observedDirection,
                        true
                );

                writeGeneRow(
                        geneWriter,
                        legacy,
                        row,
                        summary,
                        minTvd,
                        wpiThreshold,
                        analysisMode,
                        observedChange,
                        formalSupport,
                        apaChange,
                        apaDirection,
                        observedDirection,
                        apaClassification,
                        observedClassification,
                        evidenceStatus,
                        primaryP,
                        primaryQ,
                        distributionP,
                        distributionQ,
                        primaryTest,
                        distributionTest
                );

                writeSiteRows(
                        siteWriter,
                        group,
                        geneSites,
                        summary
                );
            }
        }

        ScCanonicalOutputAdapter.write(outDir, args);

        System.out.printf(
                Locale.ROOT,
                "APAURA single-cell usage harmonization complete.%n"
                        + "Gene usage: %s%n"
                        + "Site PAU: %s%n",
                geneOutput.toAbsolutePath(),
                siteOutput.toAbsolutePath()
        );
    }

    private static UsageSummary summarize(
            String group,
            List<Site> sites,
            Iterable<Sample> samples,
            Map<CountKey, Map<String, Long>> counts,
            String control,
            String treated,
            int minSampleCount
    ) {
        UsageSummary summary = new UsageSummary(sites.size());

        for (Sample sample : samples) {
            if (!sample.condition.equals(control)
                    && !sample.condition.equals(treated)) {
                continue;
            }

            Map<String, Long> sampleCounts = counts.get(
                    new CountKey(sample.id, group, sites.get(0).geneId)
            );
            if (sampleCounts == null) {
                sampleCounts = Collections.emptyMap();
            }

            long total = 0L;
            double[] pau = new double[sites.size()];
            for (int index = 0; index < sites.size(); index++) {
                long value = value(sampleCounts, sites.get(index).siteId);
                if (value < 0L) {
                    throw new IllegalArgumentException("Negative count");
                }
                total += value;
                pau[index] = value;
            }

            if (total < minSampleCount || total <= 0L) {
                continue;
            }

            for (int index = 0; index < pau.length; index++) {
                pau[index] /= total;
            }

            if (sample.condition.equals(control)) {
                addInPlace(summary.controlPau, pau);
                summary.controlSamples++;
            } else {
                addInPlace(summary.treatedPau, pau);
                summary.treatedSamples++;
            }
        }

        averageInPlace(summary.controlPau, summary.controlSamples);
        averageInPlace(summary.treatedPau, summary.treatedSamples);

        if (summary.controlSamples == 0) {
            fillNaN(summary.controlPau);
        }
        if (summary.treatedSamples == 0) {
            fillNaN(summary.treatedPau);
        }

        double[] positions = normalizedPositions(sites);
        summary.positions = positions;
        summary.controlWpi = ApaUsageMetrics.weightedPositionIndex(
                summary.controlPau,
                positions
        );
        summary.treatedWpi = ApaUsageMetrics.weightedPositionIndex(
                summary.treatedPau,
                positions
        );
        summary.deltaWpi = finiteDifference(
                summary.treatedWpi,
                summary.controlWpi
        );
        summary.apd = ApaUsageMetrics.totalVariationDistance(
                summary.controlPau,
                summary.treatedPau
        );
        summary.jsd = ApaUsageMetrics.jensenShannonDivergence(
                summary.controlPau,
                summary.treatedPau
        );
        return summary;
    }


    private static void validateAgainstLegacy(
            TsvTable legacy,
            String[] row,
            List<Site> sites,
            UsageSummary summary
    ) {
        String group = legacy.get(row, "group");
        String geneId = legacy.get(row, "gene_id");
        int legacySites = Integer.parseInt(legacy.get(row, "n_sites"));
        int legacyControlSamples = Integer.parseInt(
                legacy.get(row, "n_control_samples")
        );
        int legacyTreatedSamples = Integer.parseInt(
                legacy.get(row, "n_treated_samples")
        );

        if (legacySites != sites.size()) {
            throw new IllegalStateException(
                    "Single-cell harmonization site-count mismatch for "
                            + group + "/" + geneId
            );
        }
        if (legacyControlSamples != summary.controlSamples
                || legacyTreatedSamples != summary.treatedSamples) {
            throw new IllegalStateException(
                    "Single-cell harmonization sample-count mismatch for "
                            + group + "/" + geneId
                            + ": core=" + legacyControlSamples
                            + "/" + legacyTreatedSamples
                            + ", APAURA=" + summary.controlSamples
                            + "/" + summary.treatedSamples
            );
        }

        requireCloseWhenFinite(
                group,
                geneId,
                "control_wpi",
                parseDouble(legacy.get(row, "control_utr_index")),
                summary.controlWpi
        );
        requireCloseWhenFinite(
                group,
                geneId,
                "treated_wpi",
                parseDouble(legacy.get(row, "treated_utr_index")),
                summary.treatedWpi
        );
        requireCloseWhenFinite(
                group,
                geneId,
                "delta_wpi",
                parseDouble(legacy.get(row, "delta_utr_index")),
                summary.deltaWpi
        );
    }

    private static void requireCloseWhenFinite(
            String group,
            String geneId,
            String metric,
            double legacyValue,
            double apauraValue
    ) {
        if (!Double.isFinite(legacyValue) && !Double.isFinite(apauraValue)) {
            return;
        }
        if (!Double.isFinite(legacyValue)
                || !Double.isFinite(apauraValue)
                || Math.abs(legacyValue - apauraValue) > 5.0e-6) {
            throw new IllegalStateException(
                    "Single-cell harmonization " + metric
                            + " mismatch for " + group + "/" + geneId
                            + ": core=" + legacyValue
                            + ", APAURA=" + apauraValue
            );
        }
    }

    private static void addInPlace(double[] target, double[] values) {
        if (target.length != values.length) {
            throw new IllegalArgumentException(
                    "APA usage vectors must have the same length"
            );
        }
        for (int index = 0; index < target.length; index++) {
            target[index] += values[index];
        }
    }

    private static void averageInPlace(double[] values, int count) {
        if (count <= 0) {
            return;
        }
        for (int index = 0; index < values.length; index++) {
            values[index] /= count;
        }
    }

    private static double[] normalizedPositions(List<Site> sites) {
        double[] positions = new double[sites.size()];
        if (sites.size() <= 1) {
            return positions;
        }

        long first = sites.get(0).position;
        long last = sites.get(sites.size() - 1).position;
        double span = Math.abs((double) last - (double) first);

        for (int index = 0; index < sites.size(); index++) {
            if (span == 0.0) {
                positions[index] = 0.0;
            } else if (sites.get(0).strand == '+') {
                positions[index] =
                        (sites.get(index).position - first) / span;
            } else {
                positions[index] =
                        (first - sites.get(index).position) / span;
            }
        }
        return positions;
    }

    private static Map<String, Site> loadSites(Path path) throws IOException {
        TsvTable table = TsvTable.read(path);
        table.require(
                "site_id",
                "gene_id",
                "gene_symbol",
                "chrom",
                "position",
                "strand"
        );

        Map<String, Site> sites = new LinkedHashMap<String, Site>();
        for (String[] row : table.rows) {
            String siteId = table.get(row, "site_id");
            String strandText = table.get(row, "strand");
            if (strandText.isEmpty()) {
                throw new IllegalArgumentException(
                        "Blank strand for site " + siteId
                );
            }

            Site site = new Site(
                    siteId,
                    table.get(row, "gene_id"),
                    table.get(row, "gene_symbol"),
                    table.get(row, "chrom"),
                    Long.parseLong(table.get(row, "position")),
                    strandText.charAt(0)
            );

            if (site.strand != '+' && site.strand != '-') {
                throw new IllegalArgumentException(
                        "Invalid strand for site " + siteId
                );
            }
            if (sites.put(siteId, site) != null) {
                throw new IllegalArgumentException(
                        "Duplicate site_id: " + siteId
                );
            }
        }
        return sites;
    }

    private static Map<String, List<Site>> groupSites(
            Iterable<Site> sites
    ) {
        Map<String, List<Site>> genes =
                new LinkedHashMap<String, List<Site>>();
        for (Site site : sites) {
            List<Site> geneSites = genes.get(site.geneId);
            if (geneSites == null) {
                geneSites = new ArrayList<Site>();
                genes.put(site.geneId, geneSites);
            }
            geneSites.add(site);
        }

        for (List<Site> geneSites : genes.values()) {
            final char strand = geneSites.get(0).strand;
            for (Site site : geneSites) {
                if (site.strand != strand) {
                    throw new IllegalArgumentException(
                            "Mixed strands for gene " + site.geneId
                    );
                }
            }
            Collections.sort(geneSites, new Comparator<Site>() {
                @Override
                public int compare(Site left, Site right) {
                    int order = Long.compare(left.position, right.position);
                    return strand == '-' ? -order : order;
                }
            });
        }
        return genes;
    }

    private static Map<String, Sample> loadSamples(Path path)
            throws IOException {
        TsvTable table = TsvTable.read(path);
        table.require("sample_id", "condition");
        Map<String, Sample> samples =
                new LinkedHashMap<String, Sample>();
        for (String[] row : table.rows) {
            Sample sample = new Sample(
                    table.get(row, "sample_id"),
                    table.get(row, "condition")
            );
            if (samples.put(sample.id, sample) != null) {
                throw new IllegalArgumentException(
                        "Duplicate sample_id: " + sample.id
                );
            }
        }
        return samples;
    }

    private static void validateConditions(
            Iterable<Sample> samples,
            String control,
            String treated
    ) {
        if (control.equals(treated)) {
            throw new IllegalArgumentException(
                    "--control and --treated must differ"
            );
        }
        boolean foundControl = false;
        boolean foundTreated = false;
        for (Sample sample : samples) {
            foundControl |= sample.condition.equals(control);
            foundTreated |= sample.condition.equals(treated);
        }
        if (!foundControl) {
            throw new IllegalArgumentException(
                    "Control condition not found: " + control
            );
        }
        if (!foundTreated) {
            throw new IllegalArgumentException(
                    "Treated condition not found: " + treated
            );
        }
    }

    private static Map<CountKey, Map<String, Long>> loadCounts(
            Path countsPath,
            Map<String, Site> sites,
            Path cellsPath,
            Grouping grouping
    ) throws IOException {
        TsvTable counts = TsvTable.read(countsPath);
        if (counts.has("barcode")) {
            if (cellsPath == null) {
                throw new IllegalArgumentException(
                        "--cells is required when counts contain a barcode column"
                );
            }
            return loadCellCounts(
                    counts,
                    sites,
                    loadCellGroups(cellsPath, grouping)
            );
        }
        return loadGroupCounts(counts, sites);
    }

    private static Map<String, String> loadCellGroups(
            Path path,
            Grouping grouping
    ) throws IOException {
        TsvTable table = TsvTable.read(path);
        table.require("sample_id", "barcode");
        for (String column : grouping.columns) {
            table.require(column);
        }

        Map<String, String> groups = new LinkedHashMap<String, String>();
        for (String[] row : table.rows) {
            String key = table.get(row, "sample_id")
                    + "\t" + table.get(row, "barcode");
            String group = grouping.group(table, row);
            String previous = groups.put(key, group);
            if (previous != null && !previous.equals(group)) {
                throw new IllegalArgumentException(
                        "Conflicting metadata for " + key
                );
            }
        }
        return groups;
    }

    private static Map<CountKey, Map<String, Long>> loadCellCounts(
            TsvTable counts,
            Map<String, Site> sites,
            Map<String, String> groups
    ) {
        counts.require("sample_id", "barcode", "site_id", "count");
        Map<CountKey, Map<String, Long>> result =
                new LinkedHashMap<CountKey, Map<String, Long>>();

        for (String[] row : counts.rows) {
            String sample = counts.get(row, "sample_id");
            String barcode = counts.get(row, "barcode");
            String siteId = counts.get(row, "site_id");
            Site site = sites.get(siteId);
            if (site == null) {
                continue;
            }

            String group = groups.get(sample + "\t" + barcode);
            if (group == null) {
                continue;
            }

            long count = Long.parseLong(counts.get(row, "count"));
            if (count < 0L) {
                throw new IllegalArgumentException("Negative count");
            }
            mergeCount(
                    result,
                    new CountKey(sample, group, site.geneId),
                    siteId,
                    count
            );
        }
        return result;
    }

    private static Map<CountKey, Map<String, Long>> loadGroupCounts(
            TsvTable counts,
            Map<String, Site> sites
    ) {
        counts.require("sample_id", "group", "site_id", "count");
        Map<CountKey, Map<String, Long>> result =
                new LinkedHashMap<CountKey, Map<String, Long>>();

        for (String[] row : counts.rows) {
            String sample = counts.get(row, "sample_id");
            String group = counts.get(row, "group");
            String siteId = counts.get(row, "site_id");
            if (group.trim().isEmpty()) {
                group = "ALL";
            }
            Site site = sites.get(siteId);
            if (site == null) {
                continue;
            }

            long count = Long.parseLong(counts.get(row, "count"));
            if (count < 0L) {
                throw new IllegalArgumentException("Negative count");
            }
            mergeCount(
                    result,
                    new CountKey(sample, group, site.geneId),
                    siteId,
                    count
            );
        }
        return result;
    }

    private static void mergeCount(
            Map<CountKey, Map<String, Long>> counts,
            CountKey key,
            String siteId,
            long value
    ) {
        Map<String, Long> bySite = counts.get(key);
        if (bySite == null) {
            bySite = new LinkedHashMap<String, Long>();
            counts.put(key, bySite);
        }
        Long previous = bySite.get(siteId);
        bySite.put(siteId, previous == null ? value : previous + value);
    }

    private static long value(Map<String, Long> counts, String siteId) {
        Long value = counts.get(siteId);
        return value == null ? 0L : value;
    }

    private static void writeGeneHeader(BufferedWriter writer)
            throws IOException {
        writer.write(
                "group\tgene_id\tgene_symbol\tchrom\tstrand\tn_sites"
                        + "\tanalysis_mode\tn_control_samples\tn_treated_samples"
                        + "\tcontrol_total\ttreated_total"
                        + "\tcontrol_wpi\ttreated_wpi\tdelta_wpi"
                        + "\tapd\tusage_jsd\tmin_tvd\twpi_threshold"
                        + "\tobserved_change\tformal_support\tapa_change"
                        + "\tapa_direction\tobserved_direction"
                        + "\tapa_classification\tobserved_classification"
                        + "\tevidence_status"
                        + "\tprimary_pvalue\tprimary_qvalue"
                        + "\tdistribution_pvalue\tdistribution_qvalue"
                        + "\tprimary_test\tdistribution_test"
                        + "\tlegacy_classification"
        );
        writer.newLine();
    }

    private static void writeSiteHeader(BufferedWriter writer)
            throws IOException {
        writer.write(
                "group\tgene_id\tgene_symbol\tchrom\tstrand"
                        + "\tsite_id\tposition\tsite_rank\tnormalized_position"
                        + "\tn_control_samples\tn_treated_samples"
                        + "\tcontrol_pau\ttreated_pau\tdelta_pau"
        );
        writer.newLine();
    }

    private static void writeGeneRow(
            BufferedWriter writer,
            TsvTable legacy,
            String[] row,
            UsageSummary summary,
            double minTvd,
            double wpiThreshold,
            String analysisMode,
            boolean observedChange,
            boolean formalSupport,
            boolean apaChange,
            String apaDirection,
            String observedDirection,
            String apaClassification,
            String observedClassification,
            String evidenceStatus,
            double primaryP,
            double primaryQ,
            double distributionP,
            double distributionQ,
            String primaryTest,
            String distributionTest
    ) throws IOException {
        writer.write(String.join(
                "\t",
                legacy.get(row, "group"),
                legacy.get(row, "gene_id"),
                legacy.get(row, "gene_symbol"),
                legacy.get(row, "chrom"),
                legacy.get(row, "strand"),
                legacy.get(row, "n_sites"),
                analysisMode,
                Integer.toString(summary.controlSamples),
                Integer.toString(summary.treatedSamples),
                legacy.get(row, "control_total"),
                legacy.get(row, "treated_total"),
                format(summary.controlWpi),
                format(summary.treatedWpi),
                format(summary.deltaWpi),
                format(summary.apd),
                format(summary.jsd),
                format(minTvd),
                format(wpiThreshold),
                yesNo(observedChange),
                yesNo(formalSupport),
                yesNo(apaChange),
                apaDirection,
                observedDirection,
                apaClassification,
                observedClassification,
                evidenceStatus,
                format(primaryP),
                format(primaryQ),
                format(distributionP),
                format(distributionQ),
                primaryTest,
                distributionTest,
                legacy.get(row, "classification")
        ));
        writer.newLine();
    }

    private static void writeSiteRows(
            BufferedWriter writer,
            String group,
            List<Site> sites,
            UsageSummary summary
    ) throws IOException {
        for (int index = 0; index < sites.size(); index++) {
            Site site = sites.get(index);
            double delta = finiteDifference(
                    summary.treatedPau[index],
                    summary.controlPau[index]
            );
            writer.write(String.join(
                    "\t",
                    group,
                    site.geneId,
                    site.geneSymbol,
                    site.chrom,
                    Character.toString(site.strand),
                    site.siteId,
                    Long.toString(site.position),
                    Integer.toString(index + 1),
                    format(summary.positions[index]),
                    Integer.toString(summary.controlSamples),
                    Integer.toString(summary.treatedSamples),
                    format(summary.controlPau[index]),
                    format(summary.treatedPau[index]),
                    format(delta)
            ));
            writer.newLine();
        }
    }

    private static String evidenceStatus(
            String analysisMode,
            String primaryTest,
            boolean formalSupport
    ) {
        if ("LOW_COVERAGE".equals(primaryTest)) {
            return "low_coverage";
        }
        if (!"REPLICATED_DIFFERENTIAL".equals(analysisMode)) {
            return "descriptive_only";
        }
        return formalSupport ? "supported" : "not_significant";
    }

    private static String classification(
            boolean changed,
            String direction,
            boolean observed
    ) {
        if (!changed) {
            return observed ? "no_major_observed_shift" : "stable";
        }
        String prefix = observed ? "observed_" : "";
        if ("lengthening".equals(direction)) {
            return prefix + "lengthening";
        }
        if ("shortening".equals(direction)) {
            return prefix + "shortening";
        }
        return prefix + "redistribution";
    }

    private static boolean finiteAtMost(double value, double threshold) {
        return Double.isFinite(value) && value <= threshold;
    }

    private static double finiteDifference(double right, double left) {
        if (!Double.isFinite(right) || !Double.isFinite(left)) {
            return Double.NaN;
        }
        return right - left;
    }

    private static void fillNaN(double[] values) {
        for (int index = 0; index < values.length; index++) {
            values[index] = Double.NaN;
        }
    }

    private static String yesNo(boolean value) {
        return value ? "yes" : "no";
    }

    private static String format(double value) {
        if (!Double.isFinite(value)) {
            return "NA";
        }
        return String.format(Locale.ROOT, "%.6g", value);
    }

    private static void requireProbabilityThreshold(
            String name,
            double value
    ) {
        if (!Double.isFinite(value) || value < 0.0 || value > 1.0) {
            throw new IllegalArgumentException(
                    name + " must be between 0 and 1"
            );
        }
    }

    private static Path requiredPath(String[] args, String key) {
        return Paths.get(requiredArg(args, key));
    }

    private static Path optionalPath(String[] args, String key) {
        String value = findArg(args, key);
        return value == null ? null : Paths.get(value);
    }

    private static String requiredArg(String[] args, String key) {
        String value = findArg(args, key);
        if (value == null || value.trim().isEmpty()) {
            throw new IllegalArgumentException(
                    "Missing required option " + key
            );
        }
        return value;
    }

    private static String stringArg(
            String[] args,
            String key,
            String defaultValue
    ) {
        String value = findArg(args, key);
        return value == null ? defaultValue : value;
    }

    private static int intArg(String[] args, String key, int defaultValue) {
        String value = findArg(args, key);
        return value == null ? defaultValue : Integer.parseInt(value);
    }

    private static double doubleArg(
            String[] args,
            String key,
            double defaultValue
    ) {
        String value = findArg(args, key);
        return value == null ? defaultValue : Double.parseDouble(value);
    }

    private static boolean hasFlag(String[] args, String key) {
        for (String arg : args) {
            if (key.equals(arg)) {
                return true;
            }
        }
        return false;
    }

    private static String findArg(String[] args, String key) {
        for (int index = 0; index < args.length; index++) {
            if (!key.equals(args[index])) {
                continue;
            }
            if (index + 1 < args.length
                    && !args[index + 1].startsWith("--")) {
                return args[index + 1];
            }
            return null;
        }
        return null;
    }

    private static String groupByArg(String[] args) {
        String groupBy = findArg(args, "--group-by");
        if (groupBy != null) {
            return groupBy;
        }
        return stringArg(args, "--group-column", "cell_type");
    }

    private static double parseDouble(String value) {
        if (value == null
                || value.isEmpty()
                || "NA".equalsIgnoreCase(value)
                || "NaN".equalsIgnoreCase(value)) {
            return Double.NaN;
        }
        return Double.parseDouble(value);
    }

    private static final class Site {
        final String siteId;
        final String geneId;
        final String geneSymbol;
        final String chrom;
        final long position;
        final char strand;

        Site(
                String siteId,
                String geneId,
                String geneSymbol,
                String chrom,
                long position,
                char strand
        ) {
            this.siteId = siteId;
            this.geneId = geneId;
            this.geneSymbol = geneSymbol;
            this.chrom = chrom;
            this.position = position;
            this.strand = strand;
        }
    }

    private static final class Sample {
        final String id;
        final String condition;

        Sample(String id, String condition) {
            this.id = id;
            this.condition = condition;
        }
    }

    private static final class CountKey {
        final String sample;
        final String group;
        final String gene;

        CountKey(String sample, String group, String gene) {
            this.sample = sample;
            this.group = group;
            this.gene = gene;
        }

        @Override
        public boolean equals(Object other) {
            if (this == other) {
                return true;
            }
            if (!(other instanceof CountKey)) {
                return false;
            }
            CountKey key = (CountKey) other;
            return sample.equals(key.sample)
                    && group.equals(key.group)
                    && gene.equals(key.gene);
        }

        @Override
        public int hashCode() {
            int hash = sample.hashCode();
            hash = 31 * hash + group.hashCode();
            hash = 31 * hash + gene.hashCode();
            return hash;
        }
    }

    private static final class Grouping {
        final List<String> columns;
        final String separator;
        final boolean labelColumns;

        Grouping(String csv, String separator, boolean labelColumns) {
            this.columns = new ArrayList<String>();
            String input = csv == null || csv.trim().isEmpty()
                    ? "cell_type"
                    : csv;
            for (String part : input.split(",")) {
                String column = part.trim();
                if (!column.isEmpty()) {
                    columns.add(column);
                }
            }
            if (columns.isEmpty()) {
                throw new IllegalArgumentException(
                        "--group-by must contain at least one metadata column"
                );
            }
            this.separator = separator == null ? "::" : separator;
            this.labelColumns = labelColumns;
        }

        String group(TsvTable table, String[] row) {
            List<String> parts = new ArrayList<String>();
            for (String column : columns) {
                String value = table.get(row, column);
                if (value == null || value.trim().isEmpty()) {
                    value = "UNANNOTATED";
                }
                parts.add(labelColumns ? column + "=" + value : value);
            }
            return join(parts, separator);
        }
    }

    private static final class UsageSummary {
        final double[] controlPau;
        final double[] treatedPau;
        double[] positions;
        int controlSamples;
        int treatedSamples;
        double controlWpi = Double.NaN;
        double treatedWpi = Double.NaN;
        double deltaWpi = Double.NaN;
        double apd = Double.NaN;
        double jsd = Double.NaN;

        UsageSummary(int nSites) {
            controlPau = new double[nSites];
            treatedPau = new double[nSites];
        }
    }

    private static final class TsvTable {
        final Map<String, Integer> header;
        final List<String[]> rows;

        TsvTable(Map<String, Integer> header, List<String[]> rows) {
            this.header = header;
            this.rows = rows;
        }

        static TsvTable read(Path path) throws IOException {
            try (BufferedReader reader = Files.newBufferedReader(
                    path,
                    StandardCharsets.UTF_8
            )) {
                String headerLine = reader.readLine();
                if (headerLine == null) {
                    throw new IllegalArgumentException(
                            "Empty TSV: " + path
                    );
                }
                if (headerLine.startsWith("\uFEFF")) {
                    headerLine = headerLine.substring(1);
                }

                String[] columns = headerLine.split("\t", -1);
                Map<String, Integer> header =
                        new LinkedHashMap<String, Integer>();
                for (int index = 0; index < columns.length; index++) {
                    String column = columns[index];
                    if (header.put(column, index) != null) {
                        throw new IllegalArgumentException(
                                "Duplicate TSV column: " + column
                        );
                    }
                }

                List<String[]> rows = new ArrayList<String[]>();
                String line;
                while ((line = reader.readLine()) != null) {
                    if (!line.trim().isEmpty()) {
                        rows.add(line.split("\t", -1));
                    }
                }
                return new TsvTable(header, rows);
            }
        }

        boolean has(String column) {
            return header.containsKey(column);
        }

        void require(String... columns) {
            for (String column : columns) {
                if (!has(column)) {
                    throw new IllegalArgumentException(
                            "Missing TSV column: " + column
                    );
                }
            }
        }

        String get(String[] row, String column) {
            Integer index = header.get(column);
            if (index == null) {
                return null;
            }
            if (index >= row.length) {
                return "";
            }
            return row[index];
        }
    }

    private static String join(List<String> values, String separator) {
        StringBuilder result = new StringBuilder();
        for (int index = 0; index < values.size(); index++) {
            if (index > 0) {
                result.append(separator);
            }
            result.append(values.get(index));
        }
        return result.toString();
    }
}
