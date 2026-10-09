package org.apaura.annotation;

import org.apaura.model.Strand;
import org.apaura.model.TranscriptModel;
import org.junit.jupiter.api.Test;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

class GtfAnnotationTest {

    @Test
    void exposesTranscriptModelsFromSingleGtfParse() throws Exception {
        Path gtf = Files.createTempFile(
                "apaura-gtf-annotation-",
                ".gtf"
        );

        try {
            String text =
                    "chr1\ttest\texon\t101\t200\t.\t+\t.\t"
                            + "gene_id \"G1\"; transcript_id \"T1\"; "
                            + "gene_name \"GENE1\"; "
                            + "gene_type \"protein_coding\"; "
                            + "transcript_type \"protein_coding\";\n"
                            + "chr1\ttest\tCDS\t121\t180\t.\t+\t0\t"
                            + "gene_id \"G1\"; transcript_id \"T1\"; "
                            + "gene_name \"GENE1\"; "
                            + "gene_type \"protein_coding\"; "
                            + "transcript_type \"protein_coding\";\n"
                            + "chr1\ttest\texon\t301\t400\t.\t+\t.\t"
                            + "gene_id \"G1\"; transcript_id \"T1\"; "
                            + "gene_name \"GENE1\"; "
                            + "gene_type \"protein_coding\"; "
                            + "transcript_type \"protein_coding\";\n"
                            + "chr2\ttest\texon\t501\t600\t.\t-\t.\t"
                            + "gene_id \"G2\"; transcript_id \"T2\"; "
                            + "gene_name \"GENE2\"; "
                            + "gene_type \"lncRNA\"; "
                            + "transcript_type \"lncRNA\";\n"
                            + "chr2\ttest\texon\t701\t800\t.\t-\t.\t"
                            + "gene_id \"G2\"; transcript_id \"T2\"; "
                            + "gene_name \"GENE2\"; "
                            + "gene_type \"lncRNA\"; "
                            + "transcript_type \"lncRNA\";\n";

            Files.write(
                    gtf,
                    text.getBytes(StandardCharsets.UTF_8)
            );

            GtfAnnotation annotation =
                    GtfAnnotation.load(gtf.toString());

            List<TranscriptModel> transcripts =
                    annotation.transcripts();

            assertEquals(2, transcripts.size());

            TranscriptModel plus = transcripts.get(0);
            assertEquals("T1", plus.transcriptId);
            assertEquals("G1", plus.geneId);
            assertEquals(Strand.PLUS, plus.strand);
            assertEquals(2, plus.exons.size());
            assertTrue(plus.hasCds());
            assertEquals(300, plus.terminalExonAnchor0());
            assertEquals(399, plus.pasPos0());

            TranscriptModel minus = transcripts.get(1);
            assertEquals("T2", minus.transcriptId);
            assertEquals(Strand.MINUS, minus.strand);
            assertEquals(2, minus.exons.size());
            assertEquals(600, minus.terminalExonAnchor0());
            assertEquals(500, minus.pasPos0());

        } finally {
            Files.deleteIfExists(gtf);
        }
    }
}
