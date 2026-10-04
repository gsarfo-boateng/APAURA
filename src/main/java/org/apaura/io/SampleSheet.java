package org.apaura.io;

import org.apaura.model.Sample;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

public final class SampleSheet {
    private SampleSheet() {
    }

    public static List<Sample> read(String path) throws IOException {
        File sheetFile = new File(path).getAbsoluteFile();
        File sheetDirectory = sheetFile.getParentFile();

        try (BufferedReader reader = new BufferedReader(new FileReader(sheetFile))) {
            String headerLine = reader.readLine();
            if (headerLine == null) {
                throw new IOException("Empty sample sheet");
            }

            ColumnIndices columns = findColumns(headerLine);
            List<Sample> samples = new ArrayList<Sample>();
            String line;

            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty() || line.startsWith("#")) {
                    continue;
                }

                String[] fields = line.split("\t");
                if (fields.length <= columns.maximumIndex()) {
                    continue;
                }

                File dataFile = new File(fields[columns.path]);
                if (!dataFile.isAbsolute()) {
                    dataFile = new File(sheetDirectory, fields[columns.path]);
                }

                samples.add(new Sample(
                        fields[columns.sample],
                        fields[columns.condition],
                        dataFile.getPath()
                ));
            }

            if (samples.isEmpty()) {
                throw new IOException("No samples found in " + path);
            }
            return samples;
        }
    }

    private static ColumnIndices findColumns(String headerLine) throws IOException {
        String[] headers = headerLine.split("\t");
        int sample = -1;
        int condition = -1;
        int path = -1;

        for (int i = 0; i < headers.length; i++) {
            String header = headers[i].trim().toLowerCase(Locale.ROOT);
            if (header.equals("sample")) {
                sample = i;
            } else if (header.equals("condition") || header.equals("treatment")) {
                condition = i;
            } else if (header.equals("bam")
                    || header.equals("file")
                    || header.equals("path")
                    || header.equals("isoquant_dir")) {
                path = i;
            }
        }

        if (sample < 0 || condition < 0 || path < 0) {
            throw new IOException(
                    "Sample sheet requires columns: sample, condition/treatment, "
                            + "and bam/file/path/isoquant_dir"
            );
        }
        return new ColumnIndices(sample, condition, path);
    }

    private static final class ColumnIndices {
        final int sample;
        final int condition;
        final int path;

        ColumnIndices(int sample, int condition, int path) {
            this.sample = sample;
            this.condition = condition;
            this.path = path;
        }

        int maximumIndex() {
            return Math.max(sample, Math.max(condition, path));
        }
    }
}
