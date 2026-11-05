package org.apaura.output;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Map;

public final class CanonicalOutputWriter {

    public static final String SCHEMA_VERSION = "1.0";

    public static final String[] APA_COLUMNS = {
            "chrom",
            "start",
            "end",
            "strand",
            "gene_id",
            "gene_name",
            "biotype",
            "apa_architecture",
            "apa_direction",
            "changed_mechanisms",
            "primary_mechanism",
            "number_sites",
            "effect_size",
            "pvalue",
            "fdr",
    };

    public static final String[] EVENT_COLUMNS = {
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
            "apa_direction",
            "effect_size",
            "pvalue",
            "fdr",
    };

    public static final String[] SITE_COLUMNS = {
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
    };

    private CanonicalOutputWriter() {
    }

    public static void writeApa(
            File file,
            List<Map<String, String>> rows
    ) throws IOException {
        writeTsv(file, APA_COLUMNS, rows);
    }

    public static void writeEvents(
            File file,
            List<Map<String, String>> rows
    ) throws IOException {
        writeTsv(file, EVENT_COLUMNS, rows);
    }

    public static void writeSites(
            File file,
            List<Map<String, String>> rows
    ) throws IOException {
        writeTsv(file, SITE_COLUMNS, rows);
    }

    public static void writeRunInfo(
            File file,
            Map<String, String> values
    ) throws IOException {
        ensureParent(file);

        try (PrintWriter writer = new PrintWriter(
                new BufferedWriter(new FileWriter(file)))) {

            writer.println(
                    "output_schema_version=" + SCHEMA_VERSION
            );

            for (Map.Entry<String, String> entry : values.entrySet()) {
                if ("output_schema_version".equals(entry.getKey())) {
                    continue;
                }

                writer.println(
                        entry.getKey()
                                + "="
                                + clean(entry.getValue())
                );
            }
        }
    }

    public static File prepareDetailsDirectory(
            File outputDirectory
    ) throws IOException {
        File details = new File(
                outputDirectory,
                "details"
        );

        if (!details.exists() && !details.mkdirs()) {
            throw new IOException(
                    "Could not create details directory: "
                            + details
            );
        }

        return details;
    }

    private static void writeTsv(
            File file,
            String[] columns,
            List<Map<String, String>> rows
    ) throws IOException {
        ensureParent(file);

        try (PrintWriter writer = new PrintWriter(
                new BufferedWriter(new FileWriter(file)))) {

            writer.println(join(columns));

            for (Map<String, String> row : rows) {
                String[] values = new String[columns.length];

                for (int i = 0; i < columns.length; i++) {
                    String column = columns[i];

                    if (!row.containsKey(column)) {
                        throw new IllegalArgumentException(
                                "Missing required column '"
                                        + column
                                        + "' while writing "
                                        + file
                        );
                    }

                    values[i] = clean(row.get(column));
                }

                writer.println(join(values));
            }
        }
    }

    private static void ensureParent(
            File file
    ) throws IOException {
        File parent = file.getParentFile();

        if (parent != null
                && !parent.exists()
                && !parent.mkdirs()) {
            throw new IOException(
                    "Could not create output directory: "
                            + parent
            );
        }
    }

    private static String clean(
            String value
    ) {
        if (value == null || value.isEmpty()) {
            return "NA";
        }

        return value
                .replace('\t', ' ')
                .replace('\n', ' ')
                .replace('\r', ' ');
    }

    private static String join(
            String[] values
    ) {
        StringBuilder builder = new StringBuilder();

        for (int i = 0; i < values.length; i++) {
            if (i > 0) {
                builder.append('\t');
            }

            builder.append(values[i]);
        }

        return builder.toString();
    }
}
