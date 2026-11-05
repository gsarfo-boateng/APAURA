package org.apaura.analysis;

import org.apaura.model.Sample;
import org.apaura.model.SiteCluster;
import org.apaura.stats.ModeratedDirichletMultinomial;
import org.apaura.stats.MultipleTesting;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Gene-level compositional inference for protocol-aware 3'-end sequencing.
 *
 * <p>The whole PAS count vector is tested jointly with APAURA's existing
 * moderated Dirichlet-multinomial implementation. PAS-level beta-binomial
 * statistics remain available separately for localizing the sites that drive
 * a significant gene-level redistribution.</p>
 */
public final class ThreeSeqGeneDifferential {

    private ThreeSeqGeneDifferential() {
        // Utility class.
    }

    public static final class GeneResult {
        public final String geneId;
        public final int nSites;
        public int controlReplicates;
        public int caseReplicates;
        public double pvalue = Double.NaN;
        public double fdr = Double.NaN;
        public boolean converged = false;
        public String status = "not_tested";

        GeneResult(String geneId, int nSites) {
            this.geneId = geneId;
            this.nSites = nSites;
        }
    }

    public static final class Batch {
        public final Map<String, GeneResult> genes =
                new LinkedHashMap<String, GeneResult>();
        public double commonPrecision = Double.NaN;

        public GeneResult get(String geneId) {
            return genes.get(geneId);
        }
    }

    public static Batch test(
            List<Sample> samples,
            List<SiteCluster> sites,
            String control,
            String caseCondition
    ) {
        Batch output = new Batch();

        Map<String, List<SiteCluster>> byGene =
                groupByGene(sites);

        List<ModeratedDirichletMultinomial.Event> events =
                new ArrayList<ModeratedDirichletMultinomial.Event>();
        List<String> testedGeneIds = new ArrayList<String>();

        for (Map.Entry<String, List<SiteCluster>> entry : byGene.entrySet()) {
            String geneId = entry.getKey();
            List<SiteCluster> geneSites = entry.getValue();

            GeneResult result =
                    new GeneResult(geneId, geneSites.size());
            output.genes.put(geneId, result);

            if (geneSites.size() < 2) {
                result.status = "insufficient_sites";
                continue;
            }

            List<int[]> rows = new ArrayList<int[]>();
            List<Integer> labels = new ArrayList<Integer>();

            for (Sample sample : samples) {
                int label;
                if (control.equals(sample.condition)) {
                    label = 0;
                } else if (caseCondition.equals(sample.condition)) {
                    label = 1;
                } else {
                    continue;
                }

                int[] counts = new int[geneSites.size()];
                int total = 0;

                for (int siteIndex = 0;
                     siteIndex < geneSites.size();
                     siteIndex++) {
                    Integer value =
                            geneSites.get(siteIndex).counts.get(sample.name);
                    counts[siteIndex] =
                            value == null ? 0 : value;
                    total += counts[siteIndex];
                }

                if (total <= 0) {
                    continue;
                }

                rows.add(counts);
                labels.add(label);

                if (label == 0) {
                    result.controlReplicates++;
                } else {
                    result.caseReplicates++;
                }
            }

            if (result.controlReplicates < 2
                    || result.caseReplicates < 2) {
                result.status = "insufficient_replicate_support";
                continue;
            }

            int[][] countMatrix =
                    new int[rows.size()][geneSites.size()];
            int[] groups = new int[rows.size()];

            for (int row = 0; row < rows.size(); row++) {
                countMatrix[row] = rows.get(row);
                groups[row] = labels.get(row);
            }

            events.add(
                    new ModeratedDirichletMultinomial.Event(
                            geneId,
                            countMatrix,
                            groups
                    )
            );
            testedGeneIds.add(geneId);
            result.status = "tested";
        }

        if (events.isEmpty()) {
            return output;
        }

        ModeratedDirichletMultinomial.BatchResult batch =
                ModeratedDirichletMultinomial.test(events);

        output.commonPrecision = batch.commonPrecision;

        List<GeneResult> finite =
                new ArrayList<GeneResult>();

        for (String geneId : testedGeneIds) {
            GeneResult result = output.genes.get(geneId);
            ModeratedDirichletMultinomial.Result dm =
                    batch.get(geneId);

            if (dm == null) {
                result.status = "dm_result_missing";
                continue;
            }

            result.converged = dm.converged;
            result.pvalue = dm.pvalue;

            if (!dm.converged || !Double.isFinite(dm.pvalue)) {
                result.status = "dm_not_converged";
                continue;
            }

            result.status = "ok";
            finite.add(result);
        }

        if (!finite.isEmpty()) {
            double[] pvalues = new double[finite.size()];
            for (int index = 0; index < finite.size(); index++) {
                pvalues[index] = finite.get(index).pvalue;
            }

            double[] adjusted = MultipleTesting.bh(pvalues);
            for (int index = 0; index < finite.size(); index++) {
                finite.get(index).fdr = adjusted[index];
            }
        }

        return output;
    }

    /**
     * Conservative run-3seq gate:
     * a gene cannot remain an APA change unless the gene-level DM test
     * supports a compositional difference. Site-level statistics are retained
     * for mechanism/PAS localization.
     */
    public static void applyGate(
            List<DifferentialApa.GeneResult> genes,
            Batch dm,
            double fdrThreshold
    ) {
        for (DifferentialApa.GeneResult gene : genes) {
            GeneResult result = dm.get(gene.geneId);

            gene.compositionPvalue = result == null
                    ? Double.NaN
                    : result.pvalue;
            gene.compositionFdr = result == null
                    ? Double.NaN
                    : result.fdr;
            gene.compositionMethod =
                    result != null && result.converged
                            && Double.isFinite(result.pvalue)
                            ? "moderated-dirichlet-multinomial-LRT"
                            : "not_tested";

            if (!gene.apaChange) {
                continue;
            }

            boolean supported =
                    result != null
                            && result.converged
                            && Double.isFinite(result.fdr)
                            && result.fdr <= fdrThreshold;

            if (supported) {
                continue;
            }

            gene.apaChange = false;
            gene.changedMechanisms = "none";
            gene.mechanismDirections = "none";
            gene.apaEvent = gene.apaType + ":" + gene.apaDirection;
        }
    }

    public static void write(File file, Batch batch) throws IOException {
        File parent = file.getParentFile();
        if (parent != null
                && !parent.exists()
                && !parent.mkdirs()) {
            throw new IOException(
                    "Cannot create output directory " + parent
            );
        }

        try (PrintWriter writer = new PrintWriter(
                new BufferedWriter(new FileWriter(file)))) {
            writer.println(
                    "gene_id\tn_sites\tn_control_evaluable\tn_case_evaluable\t"
                            + "pvalue\tfdr\tconverged\tcommon_precision\tstatus"
            );

            for (GeneResult result : batch.genes.values()) {
                writer.println(
                        result.geneId
                                + "\t" + result.nSites
                                + "\t" + result.controlReplicates
                                + "\t" + result.caseReplicates
                                + "\t" + fmt(result.pvalue)
                                + "\t" + fmt(result.fdr)
                                + "\t" + result.converged
                                + "\t" + fmt(batch.commonPrecision)
                                + "\t" + result.status
                );
            }
        }
    }

    private static Map<String, List<SiteCluster>> groupByGene(
            List<SiteCluster> sites
    ) {
        Map<String, List<SiteCluster>> grouped =
                new LinkedHashMap<String, List<SiteCluster>>();

        for (SiteCluster site : sites) {
            if (site.geneId == null
                    || ".".equals(site.geneId)
                    || site.geneId.isEmpty()) {
                continue;
            }

            List<SiteCluster> geneSites =
                    grouped.get(site.geneId);

            if (geneSites == null) {
                geneSites = new ArrayList<SiteCluster>();
                grouped.put(site.geneId, geneSites);
            }

            geneSites.add(site);
        }

        return grouped;
    }

    private static String fmt(double value) {
        return Double.isFinite(value)
                ? String.format(Locale.ROOT, "%.10g", value)
                : "NA";
    }
}
