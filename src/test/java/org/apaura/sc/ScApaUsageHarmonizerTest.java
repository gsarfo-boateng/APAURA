package org.apaura.sc;

import org.junit.jupiter.api.Test;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

class ScApaUsageHarmonizerTest {

    @Test
    void detectsMultiPasRedistributionWhenDirectionalShiftCancels()
            throws Exception {
        Path dir = Files.createTempDirectory("apaura-sc-redistribution-");
        Path out = Files.createDirectories(dir.resolve("out"));

        write(dir.resolve("sites.tsv"),
                "site_id\tgene_id\tgene_symbol\tchrom\tposition\tstrand\n"
                        + "P1\tG1\tG1\tchr1\t100\t+\n"
                        + "P2\tG1\tG1\tchr1\t200\t+\n"
                        + "P3\tG1\tG1\tchr1\t300\t+\n"
                        + "P4\tG1\tG1\tchr1\t400\t+\n");
        write(dir.resolve("samples.tsv"),
                "sample_id\tcondition\n"
                        + "C1\tcontrol\nC2\tcontrol\n"
                        + "T1\ttreated\nT2\ttreated\n");
        write(dir.resolve("counts.tsv"),
                "sample_id\tgroup\tsite_id\tcount\n"
                        + sample("C1", 25, 25, 25, 25)
                        + sample("C2", 25, 25, 25, 25)
                        + sample("T1", 19, 31, 19, 31)
                        + sample("T2", 19, 31, 19, 31));
        writeLegacy(out.resolve("all_results.tsv"),
                "ALL", "G1", 4,
                "REPLICATED_DIFFERENTIAL",
                "0.50", "0.54", "0.04",
                "0.50", "0.50", "0.01", "0.01",
                "WELCH_T_SAMPLE_UTR_INDEX",
                "SAMPLE_LABEL_PERMUTATION_JS",
                "STABLE");

        ScApaUsageHarmonizer.run(new String[]{
                "compare",
                "--sites", dir.resolve("sites.tsv").toString(),
                "--samples", dir.resolve("samples.tsv").toString(),
                "--counts", dir.resolve("counts.tsv").toString(),
                "--control", "control",
                "--treated", "treated",
                "--min-tvd", "0.10",
                "--wpi-threshold", "0.10",
                "--out", out.toString()
        });

        Table row = firstDataRow(out.resolve("apaura_usage.tsv"));
        assertEquals(0.12, row.d("apd"), 1.0e-9);
        assertEquals(0.04, row.d("delta_wpi"), 1.0e-9);
        assertEquals("yes", row.s("observed_change"));
        assertEquals("yes", row.s("formal_support"));
        assertEquals("yes", row.s("apa_change"));
        assertEquals("no_directional_shift", row.s("apa_direction"));
        assertEquals("redistribution", row.s("apa_classification"));
    }

    @Test
    void averagesBiologicalSamplePauInsteadOfPoolingDepth()
            throws Exception {
        Path dir = Files.createTempDirectory("apaura-sc-sample-mean-");
        Path out = Files.createDirectories(dir.resolve("out"));

        write(dir.resolve("sites.tsv"),
                "site_id\tgene_id\tgene_symbol\tchrom\tposition\tstrand\n"
                        + "P1\tG1\tG1\tchr1\t100\t+\n"
                        + "P2\tG1\tG1\tchr1\t200\t+\n");
        write(dir.resolve("samples.tsv"),
                "sample_id\tcondition\n"
                        + "C1\tcontrol\nC2\tcontrol\n"
                        + "T1\ttreated\nT2\ttreated\n");
        write(dir.resolve("counts.tsv"),
                "sample_id\tgroup\tsite_id\tcount\n"
                        + "C1\tALL\tP1\t90\nC1\tALL\tP2\t10\n"
                        + "C2\tALL\tP1\t1\nC2\tALL\tP2\t9\n"
                        + "T1\tALL\tP1\t50\nT1\tALL\tP2\t50\n"
                        + "T2\tALL\tP1\t5\nT2\tALL\tP2\t5\n");
        writeLegacy(out.resolve("all_results.tsv"),
                "ALL", "G1", 2,
                "REPLICATED_DIFFERENTIAL",
                "0.50", "0.50", "0.00",
                "1.0", "1.0", "1.0", "1.0",
                "WELCH_T_SAMPLE_UTR_INDEX",
                "FISHER_POOLED_DESCRIPTIVE",
                "STABLE");

        ScApaUsageHarmonizer.run(new String[]{
                "compare",
                "--sites", dir.resolve("sites.tsv").toString(),
                "--samples", dir.resolve("samples.tsv").toString(),
                "--counts", dir.resolve("counts.tsv").toString(),
                "--control", "control",
                "--treated", "treated",
                "--out", out.toString()
        });

        List<String> lines = Files.readAllLines(
                out.resolve("apaura_site_usage.tsv"),
                StandardCharsets.UTF_8
        );
        String[] header = lines.get(0).split("\\t", -1);
        String[] p1 = lines.get(1).split("\\t", -1);
        int controlPau = index(header, "control_pau");
        int treatedPau = index(header, "treated_pau");

        // C1 PAU(P1)=0.9 and C2 PAU(P1)=0.1; the biological-sample mean is 0.5.
        assertEquals(0.5, Double.parseDouble(p1[controlPau]), 1.0e-9);
        assertEquals(0.5, Double.parseDouble(p1[treatedPau]), 1.0e-9);
    }

    @Test
    void cellRowsAreAggregatedWithinBiologicalSampleAndGroup()
            throws Exception {
        Path dir = Files.createTempDirectory("apaura-sc-cell-group-");
        Path out = Files.createDirectories(dir.resolve("out"));

        write(dir.resolve("sites.tsv"),
                "site_id\tgene_id\tgene_symbol\tchrom\tposition\tstrand\n"
                        + "P1\tG1\tG1\tchr1\t100\t+\n"
                        + "P2\tG1\tG1\tchr1\t200\t+\n");
        write(dir.resolve("samples.tsv"),
                "sample_id\tcondition\n"
                        + "C1\tcontrol\nC2\tcontrol\n"
                        + "T1\ttreated\nT2\ttreated\n");
        write(dir.resolve("cells.tsv"),
                "sample_id\tbarcode\tcell_type\tcell_subtype\n"
                        + "C1\tC1A\tT\tCD4\nC1\tC1B\tT\tCD4\n"
                        + "C2\tC2A\tT\tCD4\nC2\tC2B\tT\tCD4\n"
                        + "T1\tT1A\tT\tCD4\nT1\tT1B\tT\tCD4\n"
                        + "T2\tT2A\tT\tCD4\nT2\tT2B\tT\tCD4\n");
        write(dir.resolve("counts.tsv"),
                "sample_id\tbarcode\tsite_id\tcount\n"
                        + cell("C1", "C1A", 9, 1)
                        + cell("C1", "C1B", 9, 1)
                        + cell("C2", "C2A", 8, 2)
                        + cell("C2", "C2B", 8, 2)
                        + cell("T1", "T1A", 1, 9)
                        + cell("T1", "T1B", 1, 9)
                        + cell("T2", "T2A", 2, 8)
                        + cell("T2", "T2B", 2, 8));
        writeLegacy(out.resolve("all_results.tsv"),
                "T::CD4", "G1", 2,
                "REPLICATED_DIFFERENTIAL",
                "0.15", "0.85", "0.70",
                "0.001", "0.001", "0.001", "0.001",
                "WELCH_T_SAMPLE_UTR_INDEX",
                "FISHER_POOLED_DESCRIPTIVE",
                "LENGTHENING");

        ScApaUsageHarmonizer.run(new String[]{
                "compare",
                "--sites", dir.resolve("sites.tsv").toString(),
                "--samples", dir.resolve("samples.tsv").toString(),
                "--counts", dir.resolve("counts.tsv").toString(),
                "--cells", dir.resolve("cells.tsv").toString(),
                "--group-by", "cell_type,cell_subtype",
                "--control", "control",
                "--treated", "treated",
                "--out", out.toString()
        });

        Table row = firstDataRow(out.resolve("apaura_usage.tsv"));
        assertEquals("T::CD4", row.s("group"));
        assertEquals("2", row.s("n_control_samples"));
        assertEquals("2", row.s("n_treated_samples"));
        assertTrue(row.d("delta_wpi") > 0.6);
        assertEquals("lengthening", row.s("apa_direction"));
    }

    private static String sample(String id, int a, int b, int c, int d) {
        return id + "\tALL\tP1\t" + a + "\n"
                + id + "\tALL\tP2\t" + b + "\n"
                + id + "\tALL\tP3\t" + c + "\n"
                + id + "\tALL\tP4\t" + d + "\n";
    }

    private static String cell(String sample, String barcode, int p1, int p2) {
        return sample + "\t" + barcode + "\tP1\t" + p1 + "\n"
                + sample + "\t" + barcode + "\tP2\t" + p2 + "\n";
    }

    private static void writeLegacy(
            Path path,
            String group,
            String gene,
            int nSites,
            String mode,
            String controlWpi,
            String treatedWpi,
            String deltaWpi,
            String p,
            String q,
            String dp,
            String dq,
            String test,
            String distributionTest,
            String classification
    ) throws Exception {
        write(path,
                "group\tgene_id\tgene_symbol\tchrom\tstrand\tn_sites"
                        + "\tanalysis_mode\tn_control_samples\tn_treated_samples"
                        + "\tcontrol_total\ttreated_total"
                        + "\tcontrol_utr_index\ttreated_utr_index\tdelta_utr_index"
                        + "\tpvalue\tqvalue\tdistribution_pvalue\tdistribution_qvalue"
                        + "\ttest\tdistribution_test\tclassification\n"
                        + group + "\t" + gene + "\t" + gene
                        + "\tchr1\t+\t" + nSites
                        + "\t" + mode + "\t2\t2\t200\t200"
                        + "\t" + controlWpi + "\t" + treatedWpi + "\t" + deltaWpi
                        + "\t" + p + "\t" + q + "\t" + dp + "\t" + dq
                        + "\t" + test + "\t" + distributionTest
                        + "\t" + classification + "\n");
    }

    private static Table firstDataRow(Path path) throws Exception {
        List<String> lines = Files.readAllLines(path, StandardCharsets.UTF_8);
        return new Table(
                lines.get(0).split("\\t", -1),
                lines.get(1).split("\\t", -1)
        );
    }

    private static int index(String[] header, String name) {
        for (int i = 0; i < header.length; i++) {
            if (name.equals(header[i])) {
                return i;
            }
        }
        throw new AssertionError("Missing column: " + name);
    }

    private static void write(Path path, String text) throws Exception {
        Files.write(path, text.getBytes(StandardCharsets.UTF_8));
    }

    private static final class Table {
        private final String[] header;
        private final String[] row;

        Table(String[] header, String[] row) {
            this.header = header;
            this.row = row;
        }

        String s(String column) {
            return row[index(header, column)];
        }

        double d(String column) {
            return Double.parseDouble(s(column));
        }
    }
}
