package org.apaura.output;

import org.apaura.analysis.ArchitectureCatalog;
import org.apaura.analysis.ShortReadAle;
import org.apaura.analysis.ShortReadDifferential;
import org.apaura.analysis.ShortReadIntegratedApa;
import org.apaura.analysis.ShortReadIpa;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Presentation-only adapter for conventional short-read APAURA output.
 *
 * This class does not perform statistical testing or alter APA classification.
 * It maps already-computed short-read results into the canonical APAURA
 * one-row-per-gene apa.tsv schema.
 */
public final class ShortReadCanonicalAdapter {

    private ShortReadCanonicalAdapter() {
    }

    public static void writeApa(
            File file,
            List<ShortReadIntegratedApa.Result> results,
            List<ShortReadDifferential.GroupResult> tutrGroups,
            List<ShortReadAle.GeneResult> aleGenes,
            List<ShortReadIpa.SiteResult> ipaSites
    ) throws IOException {

        Map<String, ShortReadDifferential.GroupResult> tutrById =
                new LinkedHashMap<String, ShortReadDifferential.GroupResult>();

        for (ShortReadDifferential.GroupResult group : tutrGroups) {
            if (group != null && group.group != null) {
                tutrById.put(group.group.id, group);
            }
        }

        Map<String, ShortReadAle.GeneResult> aleByGene =
                new LinkedHashMap<String, ShortReadAle.GeneResult>();

        for (ShortReadAle.GeneResult gene : aleGenes) {
            if (gene != null && gene.group != null) {
                aleByGene.put(gene.group.geneId, gene);
            }
        }

        Map<String, List<ShortReadIpa.SiteResult>> ipaByGene =
                new LinkedHashMap<String, List<ShortReadIpa.SiteResult>>();

        for (ShortReadIpa.SiteResult site : ipaSites) {
            if (site == null || site.group == null) {
                continue;
            }

            List<ShortReadIpa.SiteResult> values =
                    ipaByGene.get(site.group.geneId);

            if (values == null) {
                values = new ArrayList<ShortReadIpa.SiteResult>();
                ipaByGene.put(site.group.geneId, values);
            }

            values.add(site);
        }

        List<Map<String, String>> rows =
                new ArrayList<Map<String, String>>();

        for (ShortReadIntegratedApa.Result result : results) {

            ArchitectureCatalog.Entry a = result.architecture;
            String primary = primaryMechanism(result);

            Map<String, String> row =
                    new LinkedHashMap<String, String>();

            row.put("chrom", a.chrom);
            row.put("start", Integer.toString(a.geneStart1based));
            row.put("end", Integer.toString(a.geneEnd1based));
            row.put("strand", Character.toString(a.strand.symbol()));

            row.put("gene_id", a.geneId);
            row.put("gene_name", a.geneName);
            row.put("biotype", a.biotype);

            row.put("apa_architecture", a.architectureClass);
            row.put("apa_change", result.apaChange ? "yes" : "no");
            row.put("apa_direction", userDirection(result.apaDirection));
            row.put("changed_mechanisms", result.changedMechanisms);
            row.put("primary_mechanism", primary);

            if ("none".equals(primary)) {
                row.put("number_sites", "NA");
                row.put("effect_size", "NA");
                row.put("pvalue", "NA");
                row.put("fdr", "NA");
            } else {
                row.put(
                        "number_sites",
                        numberSites(
                                result,
                                primary,
                                tutrById,
                                aleByGene,
                                ipaByGene
                        )
                );

                row.put(
                        "effect_size",
                        fmt(mechanismEffect(result, primary))
                );

                row.put(
                        "pvalue",
                        fmt(mechanismPvalue(result, primary))
                );

                row.put(
                        "fdr",
                        fmt(mechanismFdr(result, primary))
                );
            }

            row.put("status", result.status);

            rows.add(row);
        }

        CanonicalOutputWriter.writeApa(file, rows);
    }

    private static String primaryMechanism(
            ShortReadIntegratedApa.Result result
    ) {
        if (!result.apaChange
                || result.changedMechanisms == null
                || result.changedMechanisms.trim().isEmpty()
                || "none".equals(result.changedMechanisms)) {
            return "none";
        }

        String[] mechanisms =
                result.changedMechanisms.split("\\+");

        String best = mechanisms[0];

        for (int i = 1; i < mechanisms.length; i++) {
            String candidate = mechanisms[i];

            int fdrComparison = compareFinite(
                    mechanismFdr(result, candidate),
                    mechanismFdr(result, best)
            );

            if (fdrComparison < 0) {
                best = candidate;
                continue;
            }

            if (fdrComparison == 0) {
                int pComparison = compareFinite(
                        mechanismPvalue(result, candidate),
                        mechanismPvalue(result, best)
                );

                if (pComparison < 0) {
                    best = candidate;
                }
            }
        }

        return best;
    }

    private static String numberSites(
            ShortReadIntegratedApa.Result result,
            String mechanism,
            Map<String, ShortReadDifferential.GroupResult> tutrById,
            Map<String, ShortReadAle.GeneResult> aleByGene,
            Map<String, List<ShortReadIpa.SiteResult>> ipaByGene
    ) {

        if ("TUTR".equals(mechanism)) {
            ShortReadDifferential.GroupResult group =
                    tutrById.get(result.tutrRepresentativeGroup);

            if (group == null
                    || group.group == null
                    || group.group.sites == null
                    || group.group.sites.isEmpty()) {
                return "NA";
            }

            return Integer.toString(group.group.sites.size());
        }

        if ("ALE".equals(mechanism)) {
            ShortReadAle.GeneResult gene =
                    aleByGene.get(result.architecture.geneId);

            if (gene == null
                    || gene.group == null
                    || gene.group.branches == null) {
                return "NA";
            }

            Set<Integer> positions =
                    new LinkedHashSet<Integer>();

            for (ShortReadAle.Branch branch : gene.group.branches) {
                if (branch.annotatedPas0 != null) {
                    positions.addAll(branch.annotatedPas0);
                }
            }

            return positions.isEmpty()
                    ? "NA"
                    : Integer.toString(positions.size());
        }

        if ("IPA".equals(mechanism)) {
            List<ShortReadIpa.SiteResult> sites =
                    ipaByGene.get(result.architecture.geneId);

            if (sites == null || sites.isEmpty()) {
                return "NA";
            }

            Set<Integer> positions =
                    new LinkedHashSet<Integer>();

            for (ShortReadIpa.SiteResult site : sites) {
                if (site.group.annotatedPas0 != null) {
                    positions.addAll(site.group.annotatedPas0);
                }
            }

            return positions.isEmpty()
                    ? "NA"
                    : Integer.toString(positions.size());
        }

        return "NA";
    }

    private static double mechanismEffect(
            ShortReadIntegratedApa.Result result,
            String mechanism
    ) {
        if ("TUTR".equals(mechanism)) {
            return result.tutrShiftMagnitude;
        }

        if ("ALE".equals(mechanism)) {
            return result.aleSwitchIndex;
        }

        if ("IPA".equals(mechanism)) {
            return Double.isFinite(result.ipaDeltaUsage)
                    ? Math.abs(result.ipaDeltaUsage)
                    : Double.NaN;
        }

        return Double.NaN;
    }

    private static double mechanismPvalue(
            ShortReadIntegratedApa.Result result,
            String mechanism
    ) {
        if ("TUTR".equals(mechanism)) {
            return result.tutrMinPvalue;
        }

        if ("ALE".equals(mechanism)) {
            return result.aleMinPvalue;
        }

        if ("IPA".equals(mechanism)) {
            return result.ipaMinPvalue;
        }

        return Double.NaN;
    }

    private static double mechanismFdr(
            ShortReadIntegratedApa.Result result,
            String mechanism
    ) {
        if ("TUTR".equals(mechanism)) {
            return result.tutrMinFdr;
        }

        if ("ALE".equals(mechanism)) {
            return result.aleMinFdr;
        }

        if ("IPA".equals(mechanism)) {
            return result.ipaMinFdr;
        }

        return Double.NaN;
    }

    private static int compareFinite(
            double left,
            double right
    ) {
        boolean leftFinite = Double.isFinite(left);
        boolean rightFinite = Double.isFinite(right);

        if (leftFinite && rightFinite) {
            return Double.compare(left, right);
        }

        if (leftFinite) {
            return -1;
        }

        if (rightFinite) {
            return 1;
        }

        return 0;
    }

    private static String userDirection(String direction) {
        if (direction == null) {
            return "not_applicable";
        }

        if ("redistribution".equals(direction)
                || "switch".equals(direction)) {
            return "redistribution";
        }

        return direction;
    }

    private static String fmt(double value) {
        return Double.isFinite(value)
                ? String.format(java.util.Locale.ROOT, "%.10g", value)
                : "NA";
    }
}
