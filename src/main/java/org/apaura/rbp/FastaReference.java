package org.apaura.rbp;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/** FASTA accessor that uses a .fai index when available and otherwise loads sequences in memory. */
public final class FastaReference implements Closeable {
    private static final class Entry {
        final long length;
        final long offset;
        final int lineBases;
        final int lineWidth;

        Entry(long length, long offset, int lineBases, int lineWidth) {
            this.length = length;
            this.offset = offset;
            this.lineBases = lineBases;
            this.lineWidth = lineWidth;
        }
    }

    private final Map<String, Entry> index = new HashMap<String, Entry>();
    private RandomAccessFile randomAccessFile;
    private Map<String, String> inMemorySequences;

    public FastaReference(String path) throws IOException {
        File fasta = new File(path);
        File indexFile = new File(path + ".fai");

        if (indexFile.exists() && !path.endsWith(".gz")) {
            loadIndex(indexFile);
            randomAccessFile = new RandomAccessFile(fasta, "r");
        } else {
            loadIntoMemory(path);
        }
    }

    private void loadIndex(File indexFile) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(indexFile))) {
            String line;
            while ((line = reader.readLine()) != null) {
                String[] fields = line.split("\t");
                if (fields.length < 5) {
                    continue;
                }
                index.put(
                        fields[0],
                        new Entry(
                                Long.parseLong(fields[1]),
                                Long.parseLong(fields[2]),
                                Integer.parseInt(fields[3]),
                                Integer.parseInt(fields[4])
                        )
                );
            }
        }
    }

    private void loadIntoMemory(String path) throws IOException {
        inMemorySequences = new HashMap<String, String>();
        InputStream input = new FileInputStream(path);
        if (path.endsWith(".gz")) {
            input = new GZIPInputStream(input);
        }

        try (BufferedReader reader = new BufferedReader(new InputStreamReader(input))) {
            String sequenceName = null;
            StringBuilder sequence = null;
            String line;

            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    storeSequence(sequenceName, sequence);
                    sequenceName = line.substring(1).split("\\s+", 2)[0];
                    sequence = new StringBuilder();
                } else if (sequence != null) {
                    sequence.append(line.trim());
                }
            }
            storeSequence(sequenceName, sequence);
        }
    }

    private void storeSequence(String name, StringBuilder sequence) {
        if (name == null || sequence == null) {
            return;
        }
        inMemorySequences.put(
                name,
                sequence.toString().toUpperCase(Locale.ROOT)
        );
    }

    public String window(
            String chrom,
            int position0,
            char strand,
            int upstream,
            int downstream
    ) throws IOException {
        int requestedStart;
        int requestedEnd;

        if (strand == '+') {
            requestedStart = position0 - upstream;
            requestedEnd = position0 + downstream + 1;
        } else {
            requestedStart = position0 - downstream;
            requestedEnd = position0 + upstream + 1;
        }

        int leftPadding = Math.max(0, -requestedStart);
        int sequenceStart = Math.max(0, requestedStart);
        long chromosomeLength = length(chrom);
        int rightPadding = (int) Math.max(0, requestedEnd - chromosomeLength);
        int sequenceEnd = (int) Math.min(chromosomeLength, requestedEnd);

        String core = sequenceEnd > sequenceStart
                ? subsequence(chrom, sequenceStart, sequenceEnd)
                : "";

        StringBuilder result = new StringBuilder();
        appendPadding(result, leftPadding);
        result.append(core);
        appendPadding(result, rightPadding);

        String sequence = result.toString().toUpperCase(Locale.ROOT);
        return strand == '-' ? reverseComplement(sequence) : sequence;
    }

    public long length(String chrom) {
        if (inMemorySequences != null) {
            String sequence = inMemorySequences.get(chrom);
            if (sequence == null) {
                throw new IllegalArgumentException(
                        "Chromosome not found in FASTA: " + chrom
                );
            }
            return sequence.length();
        }

        Entry entry = index.get(chrom);
        if (entry == null) {
            throw new IllegalArgumentException(
                    "Chromosome not found in FASTA index: " + chrom
            );
        }
        return entry.length;
    }

    private String subsequence(String chrom, int start, int end) throws IOException {
        if (inMemorySequences != null) {
            return inMemorySequences.get(chrom).substring(start, end);
        }

        Entry entry = index.get(chrom);
        StringBuilder sequence = new StringBuilder(end - start);
        int position = start;

        while (position < end) {
            int positionWithinLine = position % entry.lineBases;
            int basesToRead = Math.min(
                    end - position,
                    entry.lineBases - positionWithinLine
            );
            long fileOffset = entry.offset
                    + (position / entry.lineBases) * (long) entry.lineWidth
                    + positionWithinLine;

            randomAccessFile.seek(fileOffset);
            byte[] bytes = new byte[basesToRead];
            randomAccessFile.readFully(bytes);
            sequence.append(new String(bytes, StandardCharsets.US_ASCII));
            position += basesToRead;
        }
        return sequence.toString();
    }

    private static void appendPadding(StringBuilder sequence, int count) {
        for (int i = 0; i < count; i++) {
            sequence.append('N');
        }
    }

    public static String reverseComplement(String sequence) {
        StringBuilder reverseComplement = new StringBuilder(sequence.length());
        for (int i = sequence.length() - 1; i >= 0; i--) {
            reverseComplement.append(complement(sequence.charAt(i)));
        }
        return reverseComplement.toString();
    }

    private static char complement(char base) {
        switch (base) {
            case 'A':
                return 'T';
            case 'C':
                return 'G';
            case 'G':
                return 'C';
            case 'T':
            case 'U':
                return 'A';
            default:
                return 'N';
        }
    }

    @Override
    public void close() throws IOException {
        if (randomAccessFile != null) {
            randomAccessFile.close();
        }
    }
}
