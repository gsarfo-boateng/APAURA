package org.apaura.rbp;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/** Parser for MEME letter-probability matrices. */
public final class MemeParser {
    private static final Pattern WIDTH_PATTERN = Pattern.compile("\\bw\\s*=\\s*(\\d+)");

    private MemeParser() {
    }

    public static List<MemeMotif> read(String path) throws IOException {
        List<MemeMotif> motifs = new ArrayList<MemeMotif>();
        String motifId = null;
        String motifName = null;

        try (BufferedReader reader = new BufferedReader(new FileReader(path))) {
            String line;
            while ((line = reader.readLine()) != null) {
                String trimmed = line.trim();

                if (trimmed.startsWith("MOTIF ")) {
                    String[] fields = trimmed.split("\\s+");
                    motifId = fields.length > 1 ? fields[1] : "motif";
                    motifName = fields.length > 2 ? fields[2] : motifId;
                    continue;
                }

                if (trimmed.startsWith("letter-probability matrix:")) {
                    int width = parseWidth(trimmed);
                    double[][] probabilities = readMatrix(reader, motifId, width);
                    motifs.add(new MemeMotif(motifId, motifName, probabilities));
                }
            }
        }

        if (motifs.isEmpty()) {
            throw new IOException(
                    "No MEME letter-probability matrices found in " + path
            );
        }
        return motifs;
    }

    private static double[][] readMatrix(
            BufferedReader reader,
            String motifId,
            int width
    ) throws IOException {
        double[][] probabilities = new double[width][4];

        for (int row = 0; row < width; row++) {
            String line = reader.readLine();
            if (line == null) {
                throw new IOException("Unexpected EOF in motif " + motifId);
            }

            String[] fields = line.trim().split("\\s+");
            if (fields.length < 4) {
                throw new IOException(
                        "Expected 4 probabilities in motif " + motifId
                );
            }

            for (int baseIndex = 0; baseIndex < 4; baseIndex++) {
                probabilities[row][baseIndex] = Double.parseDouble(fields[baseIndex]);
            }
        }
        return probabilities;
    }

    private static int parseWidth(String matrixHeader) {
        Matcher matcher = WIDTH_PATTERN.matcher(matrixHeader);
        if (!matcher.find()) {
            throw new IllegalArgumentException(
                    "MEME matrix is missing width: " + matrixHeader
            );
        }
        return Integer.parseInt(matcher.group(1));
    }
}
