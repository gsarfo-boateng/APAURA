package org.apaura.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Small dependency-free run manifest used by downstream APAURA commands.
 *
 * File format:
 *   key<TAB>value
 *
 * Paths are written as canonical absolute paths so a result directory can be
 * plotted later from a different working directory.
 */
public final class RunMetadata {

    public static final String FILE_NAME = "apaura_run.tsv";

    private RunMetadata() {
    }

    public static void write(
            File outputDirectory,
            String version,
            String analysis,
            String sampleSheet,
            String gtf,
            String control,
            String caseCondition,
            Map<String, String> extras
    ) throws IOException {

        if (outputDirectory == null) {
            throw new IllegalArgumentException("outputDirectory is null");
        }
        if (!outputDirectory.exists() && !outputDirectory.mkdirs()) {
            throw new IOException("Cannot create output directory: " + outputDirectory);
        }

        LinkedHashMap<String, String> values = new LinkedHashMap<String, String>();
        values.put("tool", "APAURA");
        values.put("version", value(version));
        values.put("analysis", value(analysis));
        values.put("samples", canonical(sampleSheet));
        values.put("gtf", canonical(gtf));
        values.put("control", value(control));
        values.put("case", value(caseCondition));

        if (extras != null) {
            for (Map.Entry<String, String> entry : extras.entrySet()) {
                if (entry.getKey() != null && !entry.getKey().trim().isEmpty()) {
                    values.put(entry.getKey().trim(), value(entry.getValue()));
                }
            }
        }

        File file = new File(outputDirectory, FILE_NAME);
        BufferedWriter writer = new BufferedWriter(new FileWriter(file));
        try {
            writer.write("key\tvalue");
            writer.newLine();

            for (Map.Entry<String, String> entry : values.entrySet()) {
                writer.write(escape(entry.getKey()));
                writer.write('\t');
                writer.write(escape(entry.getValue()));
                writer.newLine();
            }
        } finally {
            writer.close();
        }
    }

    public static Map<String, String> read(File resultDirectory) throws IOException {
        LinkedHashMap<String, String> values = new LinkedHashMap<String, String>();
        File file = new File(resultDirectory, FILE_NAME);

        if (!file.isFile()) {
            return values;
        }

        BufferedReader reader = new BufferedReader(new FileReader(file));
        try {
            String line;
            boolean first = true;

            while ((line = reader.readLine()) != null) {
                if (first) {
                    first = false;
                    if (line.equals("key\tvalue")) {
                        continue;
                    }
                }

                if (line.trim().isEmpty() || line.startsWith("#")) {
                    continue;
                }

                int tab = line.indexOf('\t');
                if (tab < 0) {
                    continue;
                }

                String key = unescape(line.substring(0, tab));
                String val = unescape(line.substring(tab + 1));
                values.put(key, val);
            }
        } finally {
            reader.close();
        }

        return values;
    }

    private static String canonical(String path) throws IOException {
        if (path == null || path.trim().isEmpty()) {
            return "";
        }
        return new File(path).getCanonicalPath();
    }

    private static String value(String value) {
        return value == null ? "" : value;
    }

    private static String escape(String value) {
        return value(value)
                .replace("\\", "\\\\")
                .replace("\t", "\\t")
                .replace("\r", "\\r")
                .replace("\n", "\\n");
    }

    private static String unescape(String value) {
        StringBuilder out = new StringBuilder();
        boolean slash = false;

        for (int i = 0; i < value.length(); i++) {
            char c = value.charAt(i);

            if (!slash) {
                if (c == '\\') {
                    slash = true;
                } else {
                    out.append(c);
                }
                continue;
            }

            if (c == 't') {
                out.append('\t');
            } else if (c == 'r') {
                out.append('\r');
            } else if (c == 'n') {
                out.append('\n');
            } else {
                out.append(c);
            }
            slash = false;
        }

        if (slash) {
            out.append('\\');
        }
        return out.toString();
    }
}
