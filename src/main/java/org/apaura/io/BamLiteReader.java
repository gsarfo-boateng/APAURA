package org.apaura.io;

import org.apaura.model.Interval;

import java.io.EOFException;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

/**
 * Sequential BAM reader used by APAURA for streaming alignment records.
 *
 * <p>The reader intentionally extracts only the fields needed by APAURA:
 * reference, position, mapping quality, flags, read length, and aligned CIGAR
 * blocks. It does not provide random access or the full BAM optional-tag API.</p>
 */
public final class BamLiteReader implements AlignmentReader {
    private static final int MIN_RECORD_SIZE = 32;
    private static final int MAX_RECORD_SIZE = 100_000_000;

    private final BgzfInputStream input;
    private final List<String> references = new ArrayList<String>();

    public BamLiteReader(String path) throws IOException {
        input = new BgzfInputStream(new FileInputStream(path));
        readHeader(path);
    }

    private void readHeader(String path) throws IOException {
        byte[] magic = new byte[4];
        readFully(magic);
        if (magic[0] != 'B'
                || magic[1] != 'A'
                || magic[2] != 'M'
                || magic[3] != 1) {
            throw new IOException("Invalid BAM magic: " + path);
        }

        int headerTextLength = readIntLittleEndian();
        skipFully(headerTextLength);

        int referenceCount = readIntLittleEndian();
        for (int index = 0; index < referenceCount; index++) {
            int nameLength = readIntLittleEndian();
            byte[] nameBytes = new byte[nameLength];
            readFully(nameBytes);
            references.add(new String(
                    nameBytes,
                    0,
                    Math.max(0, nameLength - 1),
                    StandardCharsets.UTF_8
            ));
            readIntLittleEndian(); // Reference length is not needed downstream.
        }
    }

    @Override
    public AlignmentRecord next() throws IOException {
        Integer blockSize = readIntLittleEndianOrEof();
        if (blockSize == null) {
            return null;
        }
        if (blockSize < MIN_RECORD_SIZE || blockSize > MAX_RECORD_SIZE) {
            throw new IOException("Invalid BAM record block size: " + blockSize);
        }

        byte[] data = new byte[blockSize];
        readFully(data);
        ByteBuffer buffer = ByteBuffer.wrap(data).order(ByteOrder.LITTLE_ENDIAN);

        int referenceId = buffer.getInt();
        int position0 = buffer.getInt();
        int binMapqNameLength = buffer.getInt();
        int flagCigarCount = buffer.getInt();
        int sequenceLength = buffer.getInt();
        buffer.getInt(); // mate reference id
        buffer.getInt(); // mate position
        buffer.getInt(); // template length

        int readNameLength = binMapqNameLength & 0xFF;
        int mapq = (binMapqNameLength >>> 8) & 0xFF;
        int cigarCount = flagCigarCount & 0xFFFF;
        int flag = (flagCigarCount >>> 16) & 0xFFFF;

        byte[] readNameBytes = new byte[readNameLength];
        buffer.get(readNameBytes);
        String readName = new String(
                readNameBytes,
                0,
                Math.max(0, readNameLength - 1),
                StandardCharsets.US_ASCII
        );

        CigarResult cigar = readCigar(buffer, position0, cigarCount);
        String referenceName = referenceId >= 0 && referenceId < references.size()
                ? references.get(referenceId)
                : null;

        return new AlignmentRecord(
                readName,
                referenceName,
                position0,
                mapq,
                flag,
                Math.max(0, cigar.referenceEnd0 - position0),
                sequenceLength,
                cigar.alignedBlocks
        );
    }

    private static CigarResult readCigar(
            ByteBuffer buffer,
            int position0,
            int cigarCount
    ) {
        int referencePosition = position0;
        List<Interval> alignedBlocks = new ArrayList<Interval>();

        for (int index = 0; index < cigarCount; index++) {
            int encoded = buffer.getInt();
            int length = encoded >>> 4;
            int operation = encoded & 0x0F;

            if (operation == 0 || operation == 7 || operation == 8) {
                alignedBlocks.add(new Interval(
                        referencePosition,
                        referencePosition + length
                ));
                referencePosition += length;
            } else if (operation == 2 || operation == 3) {
                referencePosition += length;
            }
        }
        return new CigarResult(referencePosition, alignedBlocks);
    }

    private static final class CigarResult {
        final int referenceEnd0;
        final List<Interval> alignedBlocks;

        CigarResult(int referenceEnd0, List<Interval> alignedBlocks) {
            this.referenceEnd0 = referenceEnd0;
            this.alignedBlocks = alignedBlocks;
        }
    }

    private void readFully(byte[] buffer) throws IOException {
        int offset = 0;
        while (offset < buffer.length) {
            int count = input.read(buffer, offset, buffer.length - offset);
            if (count < 0) {
                throw new EOFException();
            }
            offset += count;
        }
    }

    private void skipFully(int byteCount) throws IOException {
        byte[] buffer = new byte[Math.min(8192, Math.max(1, byteCount))];
        int remaining = byteCount;
        while (remaining > 0) {
            int count = input.read(buffer, 0, Math.min(remaining, buffer.length));
            if (count < 0) {
                throw new EOFException();
            }
            remaining -= count;
        }
    }

    private int readIntLittleEndian() throws IOException {
        Integer value = readIntLittleEndianOrEof();
        if (value == null) {
            throw new EOFException();
        }
        return value;
    }

    private Integer readIntLittleEndianOrEof() throws IOException {
        int byte0 = input.read();
        if (byte0 < 0) {
            return null;
        }
        int byte1 = input.read();
        int byte2 = input.read();
        int byte3 = input.read();
        if (byte1 < 0 || byte2 < 0 || byte3 < 0) {
            throw new EOFException();
        }
        return byte0
                | (byte1 << 8)
                | (byte2 << 16)
                | (byte3 << 24);
    }

    @Override
    public void close() throws IOException {
        input.close();
    }
}
