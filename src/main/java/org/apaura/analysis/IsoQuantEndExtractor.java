package org.apaura.analysis;

import org.apaura.model.EndKey;
import org.apaura.model.Strand;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;

/** Extracts transcript-oriented 3' ends from IsoQuant corrected-read BED output. */
public final class IsoQuantEndExtractor {

    private IsoQuantEndExtractor() {
        // Utility class.
    }

    public static Map<EndKey, Integer> extract(String path, String strandMode)
            throws IOException {
        File correctedBed = findCorrectedBed(new File(path));
        if (correctedBed == null) {
            throw new FileNotFoundException(
                    "Could not find *corrected_reads.bed or *.bed.gz under " + path
            );
        }

        Map<EndKey, Integer> counts = new HashMap<EndKey, Integer>();
        Set<String> seenReadIds = new HashSet<String>();
        long rows = 0;
        long endsRetained = 0;

        try (BufferedReader reader = open(correctedBed)) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty() || line.startsWith("#")) {
                    continue;
                }

                String[] fields = line.split("\t");
                if (fields.length < 3) {
                    continue;
                }

                Integer start0 = parseInt(fields[1]);
                Integer end0 = parseInt(fields[2]);
                if (start0 == null || end0 == null) {
                    continue;
                }
                rows++;

                String readId = fields.length > 3 ? fields[3] : null;
                if (readId != null && !readId.isEmpty() && !seenReadIds.add(readId)) {
                    continue;
                }

                Strand strand = fields.length > 5
                        ? Strand.parse(fields[5])
                        : Strand.UNKNOWN;
                strand = applyStrandMode(strand, strandMode);

                if (strand == Strand.UNKNOWN) {
                    throw new IOException(
                            "IsoQuant corrected BED lacks usable strand in column 6: "
                                    + correctedBed
                                    + ". Use --strand-mode plus/minus only if biologically justified."
                    );
                }

                int position0 = strand == Strand.MINUS ? start0 : end0 - 1;
                if (position0 < 0) {
                    continue;
                }

                EndKey end = new EndKey(fields[0], position0, strand);
                Integer previous = counts.get(end);
                counts.put(end, (previous == null ? 0 : previous) + 1);
                endsRetained++;
            }
        }

        System.err.println(
                "[APAURA] IsoQuant " + correctedBed.getName()
                        + ": rows=" + rows
                        + " retained_3p_ends=" + endsRetained
                        + " unique_positions=" + counts.size()
        );
        return counts;
    }

    private static BufferedReader open(File file) throws IOException {
        InputStream input = new FileInputStream(file);
        if (file.getName().endsWith(".gz")) {
            input = new GZIPInputStream(input);
        }
        return new BufferedReader(new InputStreamReader(input, StandardCharsets.UTF_8));
    }

    private static Strand applyStrandMode(Strand strand, String strandMode) {
        if ("flip".equalsIgnoreCase(strandMode)) {
            return strand.flip();
        }
        if ("plus".equalsIgnoreCase(strandMode)) {
            return Strand.PLUS;
        }
        if ("minus".equalsIgnoreCase(strandMode)) {
            return Strand.MINUS;
        }
        return strand;
    }

    private static Integer parseInt(String value) {
        try {
            return Integer.valueOf(value);
        } catch (NumberFormatException exception) {
            return null;
        }
    }

    private static File findCorrectedBed(File path) {
        if (path.isFile()) {
            return path;
        }

        File[] children = path.listFiles();
        if (children == null) {
            return null;
        }

        Arrays.sort(children, Comparator.comparing(File::getName));

        for (File child : children) {
            String name = child.getName();
            if (child.isFile()
                    && (name.contains("corrected_reads.bed")
                    || name.endsWith(".corrected_reads.bed")
                    || name.endsWith("corrected_reads.bed.gz"))) {
                return child;
            }
        }

        for (File child : children) {
            if (!child.isDirectory()) {
                continue;
            }
            File match = findCorrectedBed(child);
            if (match != null) {
                return match;
            }
        }
        return null;
    }
}
