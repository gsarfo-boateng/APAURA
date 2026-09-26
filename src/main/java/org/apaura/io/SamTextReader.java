package org.apaura.io;

import org.apaura.model.Interval;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.zip.GZIPInputStream;

/** Sequential reader for SAM and gzip-compressed SAM files. */
public final class SamTextReader implements AlignmentReader {
    private final BufferedReader reader;

    public SamTextReader(String path) throws IOException {
        InputStream input = new FileInputStream(path);
        if (path.endsWith(".gz")) {
            input = new GZIPInputStream(input);
        }
        reader = new BufferedReader(new InputStreamReader(input, StandardCharsets.UTF_8));
    }

    @Override
    public AlignmentRecord next() throws IOException {
        String line;
        while ((line = reader.readLine()) != null) {
            if (line.isEmpty() || line.charAt(0) == '@') {
                continue;
            }

            String[] fields = line.split("\t");
            if (fields.length < 6) {
                continue;
            }

            int flag = Integer.parseInt(fields[1]);
            String referenceName = "*".equals(fields[2]) ? null : fields[2];
            int position0 = "0".equals(fields[3])
                    ? -1
                    : Integer.parseInt(fields[3]) - 1;
            int mapq = Integer.parseInt(fields[4]);
            CigarInfo cigar = parseCigar(fields[5], position0);
            int queryLength = fields.length > 9 && !"*".equals(fields[9])
                    ? fields[9].length()
                    : cigar.queryConsumed;

            return new AlignmentRecord(
                    fields[0],
                    referenceName,
                    position0,
                    mapq,
                    flag,
                    cigar.referenceConsumed,
                    queryLength,
                    cigar.blocks
            );
        }
        return null;
    }

    public static final class CigarInfo {
        public final int referenceConsumed;
        public final int queryConsumed;
        public final List<Interval> blocks;

        CigarInfo(
                int referenceConsumed,
                int queryConsumed,
                List<Interval> blocks
        ) {
            this.referenceConsumed = referenceConsumed;
            this.queryConsumed = queryConsumed;
            this.blocks = blocks;
        }
    }

    public static CigarInfo parseCigar(String cigar, int position0) {
        if (cigar == null || "*".equals(cigar) || position0 < 0) {
            return new CigarInfo(
                    0,
                    0,
                    Collections.<Interval>emptyList()
            );
        }

        int operationLength = 0;
        int referencePosition = position0;
        int queryConsumed = 0;
        List<Interval> blocks = new ArrayList<Interval>();

        for (int index = 0; index < cigar.length(); index++) {
            char operation = cigar.charAt(index);
            if (Character.isDigit(operation)) {
                operationLength = operationLength * 10 + (operation - '0');
                continue;
            }
            if (operationLength <= 0) {
                operationLength = 0;
                continue;
            }

            if (operation == 'M' || operation == '=' || operation == 'X') {
                blocks.add(new Interval(
                        referencePosition,
                        referencePosition + operationLength
                ));
                referencePosition += operationLength;
                queryConsumed += operationLength;
            } else if (operation == 'D' || operation == 'N') {
                referencePosition += operationLength;
            } else if (operation == 'I' || operation == 'S') {
                queryConsumed += operationLength;
            }
            operationLength = 0;
        }

        return new CigarInfo(
                referencePosition - position0,
                queryConsumed,
                blocks
        );
    }

    public static int cigarReferenceConsumed(String cigar) {
        return parseCigar(cigar, 0).referenceConsumed;
    }

    @Override
    public void close() throws IOException {
        reader.close();
    }
}
