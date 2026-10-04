package org.apaura.io;

import java.io.IOException;
import java.util.Locale;

public final class AlignmentReaderFactory {
    private AlignmentReaderFactory() {
    }

    public static AlignmentReader open(String path) throws IOException {
        String lowerPath = path.toLowerCase(Locale.ROOT);
        if (lowerPath.endsWith(".bam")) {
            return new BamLiteReader(path);
        }
        if (lowerPath.endsWith(".sam") || lowerPath.endsWith(".sam.gz")) {
            return new SamTextReader(path);
        }
        throw new IllegalArgumentException(
                "Unsupported alignment format (use BAM or SAM): " + path
        );
    }
}
