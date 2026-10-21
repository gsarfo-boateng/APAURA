package org.apaura.sc;

import org.apaura.output.CanonicalOutputWriter;
import org.apaura.util.Citation;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Presentation-only adapter for APAURA single-cell comparison output.
 *
 * <p>The embedded sc3UTR engine and ScApaUsageHarmonizer remain solely
 * responsible for counting, replicate handling, statistical testing and
 * APA classification. This class only maps their completed output tables
 * into the canonical APAURA result layout.</p>
 */
final class ScCanonicalOutputAdapter {

    private static final String[] APA_COLUMNS = {
            "group",
            "chrom",
            "start",
            "end",
            "strand",
            "gene_id",
            "gene_name",
            "biotype",
            "apa_architecture",
            "apa_change",
            "apa_direction",
            "changed_mechanisms",
            "primary_mechanism",
            "number_sites",
            "effect_size",
            "pvalue",
            "fdr",
            "status"
    };

    private static final String[] EVENT_COLUMNS = {
            "group",
            "chrom",
            "start",
            "end",
            "strand",
            "gene_id",
            "gene_name",
            "biotype",
            "event_id",
            "mechanism",
            "number_sites",
            "apa_change",
            "apa_direction",
            "effect_size",
            "pvalue",
            "fdr",
            "status"
    };

    private static final String[] SITE_COLUMNS = {
            "group",
            "chrom",
            "position",
            "strand",
            "gene_id",
            "gene_name",
            "event_id",
            "mechanism",
            "feature_id",
            "feature_rank",
            "control_usage",
            "case_usage",
            "delta_usage",
            "pvalue",
            "fdr",
            "status"
    };

    private ScCanonicalOutputAdapter() {
    }

    static void write(Path outDir, String[] args) throws IOException {
        Path geneInput = outDir.resolve(ScApaUsageHarmonizer.GENE_OUTPUT);
        Path siteInput = outDir.resolve(ScApaUsageHarmonizer.SITE_OUTPUT);
        Path allResults = outDir.resolve("all_results.tsv");

        if (!Files.isRegularFile(geneInput)) {
            throw new IOException(
                    "Missing APAURA single-cell gene output: " + geneInput
            );
        }
        if (!Files.isRegularFile(siteInput)) {
            throw new IOException(
                    "Missing APAURA single-cell site output: " + siteInput
            );
        }

        Table genes = Table.read(geneInput);
        genes.require(
                "group",
                "gene_id",
                "gene_symbol",
                "chrom",
                "strand",
                "n_sites",
                "analysis_mode",
                "apd",
                "apa_change",
                "apa_direction",
                "evidence_status",
                "primary_pvalue",
                "primary_qvalue",
                "distribution_pvalue",
                "distribution_qvalue",
                "primary_test",
                "distribution_test"
        );

        Table sites = Table.read(siteInput);
        sites.require(
                "group",
                "gene_id",
                "gene_symbol",
                "chrom",
                "strand",
                "site_id",
                "position",
                "site_rank",
                "control_pau",
                "treated_pau",
                "delta_pau"
        );

        Map<String, Span> spans = spans(sites);

        List<String[]> apaRows = new ArrayList<String[]>();
        List<String[]> eventRows = new ArrayList<String[]>();
        List<String[]> siteRows = new ArrayList<String[]>();

        for (String[] row : genes.rows) {
            String group = genes.get(row, "group");
            String geneId = genes.get(row, "gene_id");
            Span span = spans.get(key(group, geneId));
            if (span == null) {
                throw new IllegalStateException(
                        "Missing canonical site span for "
                                + group + "/" + geneId
                );
            }

            Evidence evidence = selectEvidence(genes, row);
            String apaChange = genes.get(row, "apa_change");
            String eventId = eventId(group, geneId);

            apaRows.add(new String[]{
                    group,
                    genes.get(row, "chrom"),
                    Long.toString(span.start),
                    Long.toString(span.end),
                    genes.get(row, "strand"),
                    geneId,
                    genes.get(row, "gene_symbol"),
                    "NA",
                    "TUTR",
                    apaChange,
                    genes.get(row, "apa_direction"),
                    "yes".equals(apaChange) ? "TUTR" : "none",
                    "TUTR",
                    genes.get(row, "n_sites"),
                    genes.get(row, "apd"),
                    evidence.pvalue,
                    evidence.fdr,
                    genes.get(row, "evidence_status")
            });

            eventRows.add(new String[]{
                    group,
                    genes.get(row, "chrom"),
                    Long.toString(span.start),
                    Long.toString(span.end),
                    genes.get(row, "strand"),
                    geneId,
                    genes.get(row, "gene_symbol"),
                    "NA",
                    eventId,
                    "TUTR",
                    genes.get(row, "n_sites"),
                    apaChange,
                    genes.get(row, "apa_direction"),
                    genes.get(row, "apd"),
                    evidence.pvalue,
                    evidence.fdr,
                    genes.get(row, "evidence_status")
            });
        }

        for (String[] row : sites.rows) {
            String group = sites.get(row, "group");
            String geneId = sites.get(row, "gene_id");

            siteRows.add(new String[]{
                    group,
                    sites.get(row, "chrom"),
                    sites.get(row, "position"),
                    sites.get(row, "strand"),
                    geneId,
                    sites.get(row, "gene_symbol"),
                    eventId(group, geneId),
                    "TUTR",
                    sites.get(row, "site_id"),
                    sites.get(row, "site_rank"),
                    sites.get(row, "control_pau"),
                    sites.get(row, "treated_pau"),
                    sites.get(row, "delta_pau"),
                    "NA",
                    "NA",
                    "usage_only"
            });
        }

        writeTsv(outDir.resolve("apa.tsv"), APA_COLUMNS, apaRows);
        writeTsv(
                outDir.resolve("apa_events.tsv"),
                EVENT_COLUMNS,
                eventRows
        );
        writeTsv(
                outDir.resolve("apa_sites.tsv"),
                SITE_COLUMNS,
                siteRows
        );

        Path details = outDir.resolve("details");
        Files.createDirectories(details);

        copyIfPresent(
                geneInput,
                details.resolve(ScApaUsageHarmonizer.GENE_OUTPUT)
        );
        copyIfPresent(
                siteInput,
                details.resolve(ScApaUsageHarmonizer.SITE_OUTPUT)
        );
        copyIfPresent(
                allResults,
                details.resolve("all_results.tsv")
        );

        Path existingRunInfo = outDir.resolve("run_info.txt");
        if (Files.isRegularFile(existingRunInfo)
                && !isCanonicalRunInfo(existingRunInfo)) {
            copyIfPresent(
                    existingRunInfo,
                    details.resolve("sc3utr_run_info.txt")
            );
        }

        LinkedHashMap<String, String> info =
                new LinkedHashMap<String, String>();

        info.put("apaura_version", Citation.VERSION);
        info.put("output_schema_variant", "single-cell-grouped");
        info.put("mode", "sc");
        info.put("analysis", "compare");
        info.put("sites", arg(args, "--sites", "NA"));
        info.put("samples", arg(args, "--samples", "NA"));
        info.put("counts", arg(args, "--counts", "NA"));
        info.put("cells", arg(args, "--cells", "NA"));
        info.put("control", arg(args, "--control", "NA"));
        info.put("case", arg(args, "--treated", "NA"));
        info.put(
                "group_by",
                firstNonNull(
                        argOrNull(args, "--group-by"),
                        argOrNull(args, "--group-column"),
                        "cell_type"
                )
        );
        info.put(
                "coordinate_system",
                "positions copied unchanged from the single-cell --sites table; "
                        + "start/end are the analyzed PAS span"
        );
        info.put("mechanism_scope", "TUTR");
        info.put(
                "min_sample_count",
                arg(args, "--min-sample-count", "5")
        );

        String legacyMinDelta = arg(args, "--min-delta", "0.10");
        info.put(
                "effect_threshold",
                arg(args, "--min-tvd", legacyMinDelta)
        );
        info.put(
                "wpi_threshold",
                arg(args, "--wpi-threshold", legacyMinDelta)
        );
        info.put("fdr_threshold", arg(args, "--fdr", "0.05"));
        info.put(
                "effect_size",
                "APD_total_variation_distance_from_existing_apaura_usage"
        );
        info.put(
                "canonical_pvalue_fdr",
                "for REPLICATED_DIFFERENTIAL rows, choose the p/q pair with "
                        + "the lower FDR between the primary test and the formal "
                        + "SAMPLE_LABEL_PERMUTATION_JS distribution test; "
                        + "otherwise NA"
        );
        info.put(
                "classification",
                "copied unchanged from ScApaUsageHarmonizer apa_change and apa_direction"
        );
        info.put(
                "site_statistics",
                "not provided by the existing single-cell harmonizer; "
                        + "apa_sites.tsv pvalue/fdr are NA"
        );
        info.put(
                "legacy_outputs",
                "preserved in details/ and retained at the top level for backward compatibility"
        );

        CanonicalOutputWriter.writeRunInfo(
                outDir.resolve("run_info.txt").toFile(),
                info
        );
    }

    private static Evidence selectEvidence(
            Table table,
            String[] row
    ) {
        if (!"REPLICATED_DIFFERENTIAL".equals(
                table.get(row, "analysis_mode"))) {
            return Evidence.na();
        }

        Evidence primary = Evidence.of(
                table.get(row, "primary_pvalue"),
                table.get(row, "primary_qvalue")
        );

        Evidence distribution = Evidence.na();
        if ("SAMPLE_LABEL_PERMUTATION_JS".equals(
                table.get(row, "distribution_test"))) {
            distribution = Evidence.of(
                    table.get(row, "distribution_pvalue"),
                    table.get(row, "distribution_qvalue")
            );
        }

        if (!primary.valid) {
            return distribution;
        }
        if (!distribution.valid) {
            return primary;
        }

        return distribution.q < primary.q
                ? distribution
                : primary;
    }

    private static Map<String, Span> spans(Table sites) {
        Map<String, Span> result =
                new LinkedHashMap<String, Span>();

        for (String[] row : sites.rows) {
            String group = sites.get(row, "group");
            String geneId = sites.get(row, "gene_id");
            long position = Long.parseLong(
                    sites.get(row, "position")
            );

            String key = key(group, geneId);
            Span span = result.get(key);
            if (span == null) {
                span = new Span(position);
                result.put(key, span);
            } else {
                span.include(position);
            }
        }

        return result;
    }

    private static String key(String group, String geneId) {
        return group + "\t" + geneId;
    }

    private static String eventId(
            String group,
            String geneId
    ) {
        return group + "|" + geneId + "|TUTR";
    }

    private static void writeTsv(
            Path path,
            String[] header,
            List<String[]> rows
    ) throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(
                path,
                StandardCharsets.UTF_8
        )) {
            writer.write(join(header));
            writer.newLine();

            for (String[] row : rows) {
                if (row.length != header.length) {
                    throw new IllegalStateException(
                            "Canonical single-cell row has "
                                    + row.length
                                    + " fields; expected "
                                    + header.length
                    );
                }

                writer.write(join(row));
                writer.newLine();
            }
        }
    }

    private static String join(String[] values) {
        StringBuilder builder = new StringBuilder();

        for (int index = 0; index < values.length; index++) {
            if (index > 0) {
                builder.append('\t');
            }
            builder.append(clean(values[index]));
        }

        return builder.toString();
    }

    private static String clean(String value) {
        if (value == null || value.isEmpty()) {
            return "NA";
        }

        return value
                .replace('\t', ' ')
                .replace('\n', ' ')
                .replace('\r', ' ');
    }

    private static void copyIfPresent(
            Path source,
            Path destination
    ) throws IOException {
        if (!Files.isRegularFile(source)) {
            return;
        }

        Files.copy(
                source,
                destination,
                StandardCopyOption.REPLACE_EXISTING
        );
    }

    private static boolean isCanonicalRunInfo(
            Path path
    ) throws IOException {
        List<String> lines = Files.readAllLines(
                path,
                StandardCharsets.UTF_8
        );

        return !lines.isEmpty()
                && lines.get(0).startsWith(
                        "output_schema_version="
                );
    }

    private static String arg(
            String[] args,
            String key,
            String defaultValue
    ) {
        String value = argOrNull(args, key);
        return value == null ? defaultValue : value;
    }

    private static String argOrNull(
            String[] args,
            String key
    ) {
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

    private static String firstNonNull(
            String first,
            String second,
            String fallback
    ) {
        if (first != null) {
            return first;
        }
        if (second != null) {
            return second;
        }
        return fallback;
    }

    private static final class Span {
        long start;
        long end;

        Span(long position) {
            start = position;
            end = position;
        }

        void include(long position) {
            start = Math.min(start, position);
            end = Math.max(end, position);
        }
    }

    private static final class Evidence {
        final boolean valid;
        final String pvalue;
        final String fdr;
        final double q;

        Evidence(
                boolean valid,
                String pvalue,
                String fdr,
                double q
        ) {
            this.valid = valid;
            this.pvalue = pvalue;
            this.fdr = fdr;
            this.q = q;
        }

        static Evidence of(String pvalue, String fdr) {
            double p = parse(pvalue);
            double q = parse(fdr);

            if (!Double.isFinite(p) || !Double.isFinite(q)) {
                return na();
            }

            return new Evidence(
                    true,
                    format(p),
                    format(q),
                    q
            );
        }

        static Evidence na() {
            return new Evidence(
                    false,
                    "NA",
                    "NA",
                    Double.NaN
            );
        }
    }

    private static double parse(String value) {
        if (value == null
                || value.isEmpty()
                || "NA".equalsIgnoreCase(value)
                || "NaN".equalsIgnoreCase(value)) {
            return Double.NaN;
        }

        return Double.parseDouble(value);
    }

    private static String format(double value) {
        if (!Double.isFinite(value)) {
            return "NA";
        }

        return String.format(
                Locale.ROOT,
                "%.10g",
                value
        );
    }

    private static final class Table {
        final Map<String, Integer> index =
                new LinkedHashMap<String, Integer>();
        final List<String[]> rows =
                new ArrayList<String[]>();

        static Table read(Path path) throws IOException {
            List<String> lines = Files.readAllLines(
                    path,
                    StandardCharsets.UTF_8
            );

            if (lines.isEmpty()) {
                throw new IOException(
                        "Empty TSV file: " + path
                );
            }

            Table table = new Table();
            String[] header = lines.get(0).split("\\t", -1);

            for (int index = 0; index < header.length; index++) {
                table.index.put(header[index], index);
            }

            for (int line = 1; line < lines.size(); line++) {
                if (lines.get(line).trim().isEmpty()) {
                    continue;
                }

                String[] row = lines.get(line).split("\\t", -1);
                if (row.length != header.length) {
                    throw new IOException(
                            "Malformed TSV row in "
                                    + path
                                    + " at line "
                                    + (line + 1)
                    );
                }

                table.rows.add(row);
            }

            return table;
        }

        void require(String... columns) {
            for (String column : columns) {
                if (!index.containsKey(column)) {
                    throw new IllegalArgumentException(
                            "Missing required column: " + column
                    );
                }
            }
        }

        String get(String[] row, String column) {
            Integer position = index.get(column);
            if (position == null) {
                throw new IllegalArgumentException(
                        "Missing column: " + column
                );
            }

            return row[position];
        }
    }
}
