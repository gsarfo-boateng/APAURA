package org.apaura.analysis;

import org.apaura.io.AlignmentReader;
import org.apaura.io.AlignmentReaderFactory;
import org.apaura.io.AlignmentRecord;
import org.apaura.model.EndKey;
import org.apaura.model.Strand;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/** Extracts transcript-oriented 3' aligned ends from SAM/BAM input. */
public final class EndExtractor {

    private EndExtractor() {
        // Utility class.
    }

    public static Map<EndKey, Integer> extract(
            String path,
            int minMapq,
            boolean keepDuplicates,
            String strandMode
    ) throws IOException {
        Map<EndKey, Integer> counts = new HashMap<EndKey, Integer>();
        long alignmentsSeen = 0;
        long endsRetained = 0;

        try (AlignmentReader reader = AlignmentReaderFactory.open(path)) {
            AlignmentRecord alignment;
            while ((alignment = reader.next()) != null) {
                alignmentsSeen++;

                if (!isUsable(alignment, minMapq, keepDuplicates)) {
                    continue;
                }

                Strand strand = applyStrandMode(
                        Strand.fromFlag(alignment.isNegative()),
                        strandMode
                );
                int position0 = strand == Strand.MINUS
                        ? alignment.pos0
                        : alignment.end0Exclusive() - 1;

                if (position0 < 0) {
                    continue;
                }

                EndKey end = new EndKey(alignment.referenceName, position0, strand);
                Integer previous = counts.get(end);
                counts.put(end, previous == null ? 1 : previous + 1);
                endsRetained++;
            }
        }

        System.err.println(
                "[APAURA] " + new File(path).getName()
                        + ": scanned=" + alignmentsSeen
                        + " retained_3p_ends=" + endsRetained
                        + " unique_positions=" + counts.size()
        );
        return counts;
    }

    private static boolean isUsable(
            AlignmentRecord alignment,
            int minMapq,
            boolean keepDuplicates
    ) {
        if (alignment.isUnmapped()
                || alignment.isSecondary()
                || alignment.isSupplementary()
                || alignment.isQcFail()) {
            return false;
        }
        if (!keepDuplicates && alignment.isDuplicate()) {
            return false;
        }
        return alignment.mapq >= minMapq && alignment.referenceConsumed > 0;
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
}
