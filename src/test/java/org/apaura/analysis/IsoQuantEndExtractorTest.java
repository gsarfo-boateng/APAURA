package org.apaura.analysis;

import org.apaura.model.EndKey;
import org.junit.jupiter.api.Test;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;

class IsoQuantEndExtractorTest {

    @Test
    void ignoresMacAppleDoubleCorrectedBedSidecar() throws Exception {
        Path dir = Files.createTempDirectory("apaura-isoquant-appledouble-");

        Files.write(
                dir.resolve("._OUT.corrected_reads.bed"),
                new byte[0]
        );

        Files.write(
                dir.resolve("OUT.corrected_reads.bed"),
                ("chr1\t100\t151\tread1\t0\t+\n")
                        .getBytes(StandardCharsets.UTF_8)
        );

        Map<EndKey, Integer> counts =
                IsoQuantEndExtractor.extract(dir.toString(), "auto");

        assertEquals(1, counts.size());

        int total = 0;
        for (Integer count : counts.values()) {
            total += count;
        }
        assertEquals(1, total);
    }
}
