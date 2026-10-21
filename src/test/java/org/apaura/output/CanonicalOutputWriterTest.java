package org.apaura.output;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.nio.file.Files;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class CanonicalOutputWriterTest {

    @TempDir
    File tempDir;

    @Test
    void writesCanonicalApaHeader() throws Exception {
        File file = new File(tempDir, "apa.tsv");

        CanonicalOutputWriter.writeApa(
                file,
                Collections.<Map<String, String>>emptyList()
        );

        List<String> lines = Files.readAllLines(
                file.toPath()
        );

        assertEquals(
                "chrom\tstart\tend\tstrand\tgene_id\tgene_name\tbiotype\t"
                        + "apa_architecture\tapa_change\tapa_direction\t"
                        + "changed_mechanisms\tprimary_mechanism\t"
                        + "number_sites\teffect_size\tpvalue\tfdr\tstatus",
                lines.get(0)
        );
    }

    @Test
    void writesCanonicalEventHeader() throws Exception {
        File file = new File(
                tempDir,
                "apa_events.tsv"
        );

        CanonicalOutputWriter.writeEvents(
                file,
                Collections.<Map<String, String>>emptyList()
        );

        List<String> lines = Files.readAllLines(
                file.toPath()
        );

        assertEquals(
                "chrom\tstart\tend\tstrand\tgene_id\tgene_name\tbiotype\t"
                        + "event_id\tmechanism\tnumber_sites\tapa_change\t"
                        + "apa_direction\teffect_size\tpvalue\tfdr\tstatus",
                lines.get(0)
        );
    }

    @Test
    void writesCanonicalSiteHeader() throws Exception {
        File file = new File(
                tempDir,
                "apa_sites.tsv"
        );

        CanonicalOutputWriter.writeSites(
                file,
                Collections.<Map<String, String>>emptyList()
        );

        List<String> lines = Files.readAllLines(
                file.toPath()
        );

        assertEquals(
                "chrom\tposition\tstrand\tgene_id\tgene_name\tevent_id\t"
                        + "mechanism\tfeature_id\tfeature_rank\t"
                        + "control_usage\tcase_usage\tdelta_usage\t"
                        + "pvalue\tfdr\tstatus",
                lines.get(0)
        );
    }

    @Test
    void writesRunInfoSchemaVersionAndValues() throws Exception {
        File file = new File(
                tempDir,
                "run_info.txt"
        );

        Map<String, String> info =
                new LinkedHashMap<String, String>();

        info.put("apaura_version", "0.10.1");
        info.put("mode", "short-read");
        info.put("control", "control");
        info.put("case", "case");

        CanonicalOutputWriter.writeRunInfo(
                file,
                info
        );

        List<String> lines = Files.readAllLines(
                file.toPath()
        );

        assertEquals(
                "output_schema_version=1.0",
                lines.get(0)
        );

        assertTrue(
                lines.contains(
                        "apaura_version=0.10.1"
                )
        );

        assertTrue(
                lines.contains(
                        "mode=short-read"
                )
        );
    }

    @Test
    void preparesDetailsDirectory() throws Exception {
        File details =
                CanonicalOutputWriter.prepareDetailsDirectory(
                        tempDir
                );

        assertTrue(details.isDirectory());

        assertEquals(
                "details",
                details.getName()
        );
    }
}
