package org.apaura.io;

import java.io.BufferedInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

/**
 * Sequential BGZF decompressor used by {@link BamLiteReader}.
 *
 * <p>This implementation is intentionally limited to streaming BGZF blocks;
 * it does not implement virtual offsets or indexed random access.</p>
 */
public final class BgzfInputStream extends InputStream {
    private static final int GZIP_MAGIC_1 = 31;
    private static final int GZIP_MAGIC_2 = 139;

    private final BufferedInputStream input;
    private byte[] uncompressedBlock = new byte[0];
    private int blockOffset;
    private boolean eof;

    public BgzfInputStream(InputStream input) {
        this.input = new BufferedInputStream(input, 1 << 16);
    }

    @Override
    public int read() throws IOException {
        while (blockOffset >= uncompressedBlock.length) {
            if (!readNextBlock()) {
                return -1;
            }
        }
        return uncompressedBlock[blockOffset++] & 0xFF;
    }

    @Override
    public int read(byte[] buffer, int offset, int length) throws IOException {
        if (length == 0) {
            return 0;
        }

        int totalRead = 0;
        int remaining = length;
        int outputOffset = offset;

        while (remaining > 0) {
            if (blockOffset >= uncompressedBlock.length && !readNextBlock()) {
                return totalRead == 0 ? -1 : totalRead;
            }

            int count = Math.min(
                    remaining,
                    uncompressedBlock.length - blockOffset
            );
            System.arraycopy(
                    uncompressedBlock,
                    blockOffset,
                    buffer,
                    outputOffset,
                    count
            );
            blockOffset += count;
            outputOffset += count;
            remaining -= count;
            totalRead += count;
        }
        return totalRead;
    }

    private boolean readNextBlock() throws IOException {
        if (eof) {
            return false;
        }

        int id1 = input.read();
        if (id1 < 0) {
            eof = true;
            return false;
        }
        int id2 = input.read();
        if (id1 != GZIP_MAGIC_1 || id2 != GZIP_MAGIC_2) {
            throw new IOException("Not a BGZF/GZIP stream");
        }

        byte[] fixedHeader = new byte[8];
        readFully(input, fixedHeader, 0, fixedHeader.length);
        int flags = fixedHeader[1] & 0xFF;
        if ((flags & 4) == 0) {
            throw new IOException("GZIP stream lacks BGZF extra field");
        }

        byte[] extraLengthBytes = new byte[2];
        readFully(input, extraLengthBytes, 0, 2);
        int extraLength = unsigned16(extraLengthBytes, 0);
        byte[] extra = new byte[extraLength];
        readFully(input, extra, 0, extraLength);

        int blockSizeMinusOne = findBgzfBlockSize(extra);
        if (blockSizeMinusOne < 0) {
            throw new IOException("Missing BGZF BC subfield");
        }

        int totalBlockSize = blockSizeMinusOne + 1;
        int headerSize = 12 + extraLength;
        int remainingSize = totalBlockSize - headerSize;
        if (remainingSize < 8) {
            throw new IOException("Invalid BGZF block size");
        }

        byte[] remainingBlock = new byte[remainingSize];
        readFully(input, remainingBlock, 0, remainingSize);
        int compressedLength = remainingSize - 8;
        int expectedUncompressedSize = (int) unsigned32(
                remainingBlock,
                compressedLength + 4
        );

        if (expectedUncompressedSize == 0) {
            uncompressedBlock = new byte[0];
            blockOffset = 0;
            eof = true;
            return false;
        }

        uncompressedBlock = inflate(
                remainingBlock,
                compressedLength,
                expectedUncompressedSize
        );
        blockOffset = 0;
        return true;
    }

    private static int findBgzfBlockSize(byte[] extra) throws IOException {
        int offset = 0;
        while (offset + 4 <= extra.length) {
            int subfieldId1 = extra[offset] & 0xFF;
            int subfieldId2 = extra[offset + 1] & 0xFF;
            int subfieldLength = unsigned16(extra, offset + 2);
            if (offset + 4 + subfieldLength > extra.length) {
                throw new IOException("Malformed BGZF extra field");
            }
            if (subfieldId1 == 'B' && subfieldId2 == 'C' && subfieldLength == 2) {
                return unsigned16(extra, offset + 4);
            }
            offset += 4 + subfieldLength;
        }
        return -1;
    }

    private static byte[] inflate(
            byte[] compressed,
            int compressedLength,
            int expectedSize
    ) throws IOException {
        Inflater inflater = new Inflater(true);
        inflater.setInput(compressed, 0, compressedLength);
        byte[] output = new byte[expectedSize];

        try {
            int actualSize = inflater.inflate(output);
            if (actualSize != expectedSize || !inflater.finished()) {
                throw new IOException(
                        "BGZF inflate size mismatch: expected "
                                + expectedSize
                                + " got "
                                + actualSize
                );
            }
            return output;
        } catch (DataFormatException exception) {
            throw new IOException("Corrupt BGZF block", exception);
        } finally {
            inflater.end();
        }
    }

    private static int unsigned16(byte[] bytes, int offset) {
        return (bytes[offset] & 0xFF)
                | ((bytes[offset + 1] & 0xFF) << 8);
    }

    private static long unsigned32(byte[] bytes, int offset) {
        return ((long) bytes[offset] & 0xFF)
                | (((long) bytes[offset + 1] & 0xFF) << 8)
                | (((long) bytes[offset + 2] & 0xFF) << 16)
                | (((long) bytes[offset + 3] & 0xFF) << 24);
    }

    private static void readFully(
            InputStream input,
            byte[] buffer,
            int offset,
            int length
    ) throws IOException {
        int currentOffset = offset;
        int remaining = length;
        while (remaining > 0) {
            int count = input.read(buffer, currentOffset, remaining);
            if (count < 0) {
                throw new EOFException("Unexpected EOF in BGZF block");
            }
            currentOffset += count;
            remaining -= count;
        }
    }

    @Override
    public void close() throws IOException {
        input.close();
    }
}
