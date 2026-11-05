package org.apaura.analysis;

import org.apaura.model.Sample;
import org.apaura.model.ShortReadGroup;
import org.apaura.model.Strand;
import org.apaura.model.TranscriptModel;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class ShortReadNovelPasCollectorTest {

    @TempDir
    Path tempDir;

    @Test
    public void discoversAndIntegratesReplicatedNovelTutrPas()
            throws Exception {

        Path s1 = tempDir.resolve("s1.sam");
        Path s2 = tempDir.resolve("s2.sam");
        Path s3 = tempDir.resolve("s3.sam");

        writeBreakpointSam(s1, 60);
        writeBreakpointSam(s2, 61);
        writeFlatSam(s3);

        List<Sample> samples = new ArrayList<Sample>();
        samples.add(new Sample(
                "S1",
                "control",
                s1.toString()
        ));
        samples.add(new Sample(
                "S2",
                "treated",
                s2.toString()
        ));
        samples.add(new Sample(
                "S3",
                "control",
                s3.toString()
        ));

        TranscriptModel transcript =
                new TranscriptModel(
                        "TX1",
                        "GENE1",
                        "GENE1",
                        "protein_coding",
                        "chr1",
                        Strand.PLUS
                );

        transcript.addExon(100, 200);
        transcript.addCds(100, 130);

        List<TranscriptModel> transcripts =
                new ArrayList<TranscriptModel>();

        transcripts.add(transcript);

        ShortReadNovelPasCollector.Result discovery =
                ShortReadNovelPasCollector.discover(
                        samples,
                        transcripts,
                        5,
                        24,
                        10,
                        false,
                        "unstranded",
                        20,
                        10.0,
                        0.50,
                        2,
                        4,
                        20
                );

        assertEquals(1, discovery.sites.size());
        assertEquals(1, discovery.acceptedByFamily.size());

        ShortReadNovelPasCollector.Site site =
                discovery.sites.get(0);

        assertEquals("accepted", site.status);
        assertEquals(2, site.supportingSamples);
        assertTrue(site.position0 >= 158);
        assertTrue(site.position0 <= 161);
        assertTrue(site.coverageDrop > 0.75);

        ShortReadGroupBuilder.BuildResult built =
                ShortReadGroupBuilder.buildWithReport(
                        transcripts,
                        5,
                        24,
                        20,
                        discovery.acceptedByFamily
                );

        assertEquals(1, built.groups.size());

        ShortReadGroup group =
                built.groups.get(0);

        assertEquals(2, group.sites.size());

        assertTrue(
                group.sites.get(0).pos0 >= 158
                        && group.sites.get(0).pos0 <= 161
        );

        assertEquals(
                199,
                group.sites.get(1).pos0
        );
    }

    private static void writeBreakpointSam(
            Path path,
            int breakpointLength
    ) throws IOException {

        try (BufferedWriter writer =
                     Files.newBufferedWriter(
                             path,
                             StandardCharsets.UTF_8
                     )) {

            writer.write("@HD\tVN:1.6\tSO:unsorted\n");
            writer.write("@SQ\tSN:chr1\tLN:1000\n");

            /*
             * 30 reads support the upstream portion only;
             * 5 reads traverse the complete terminal region.
             *
             * The detector never examines the endpoint counts directly.
             * These alignments are converted to CIGAR-aware base coverage.
             */
            for (int i = 0; i < 30; i++) {
                writeRead(
                        writer,
                        "short" + i,
                        breakpointLength
                );
            }

            for (int i = 0; i < 5; i++) {
                writeRead(
                        writer,
                        "long" + i,
                        100
                );
            }
        }
    }

    private static void writeFlatSam(
            Path path
    ) throws IOException {

        try (BufferedWriter writer =
                     Files.newBufferedWriter(
                             path,
                             StandardCharsets.UTF_8
                     )) {

            writer.write("@HD\tVN:1.6\tSO:unsorted\n");
            writer.write("@SQ\tSN:chr1\tLN:1000\n");

            for (int i = 0; i < 35; i++) {
                writeRead(
                        writer,
                        "flat" + i,
                        100
                );
            }
        }
    }

    private static void writeRead(
            BufferedWriter writer,
            String name,
            int length
    ) throws IOException {

        writer.write(name);
        writer.write("\t0\tchr1\t101\t60\t");
        writer.write(Integer.toString(length));
        writer.write("M\t*\t0\t0\t");
        writer.write(repeat('A', length));
        writer.write("\t");
        writer.write(repeat('I', length));
        writer.write("\n");
    }

    private static String repeat(
            char value,
            int length
    ) {
        StringBuilder result =
                new StringBuilder(length);

        for (int i = 0; i < length; i++) {
            result.append(value);
        }

        return result.toString();
    }
}
