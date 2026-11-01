package org.apaura.output;

import org.apaura.analysis.ArchitectureAnalysis;
import org.apaura.analysis.DifferentialApa;
import org.apaura.model.SiteCluster;
import org.apaura.util.Citation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Presentation-only adapter for the existing direct-RNA / IsoQuant APA model.
 *
 * <p>This class does not perform statistical tests or alter APA classification.
 * It maps already-computed DifferentialApa and ArchitectureAnalysis values into
 * the canonical APAURA output schema.</p>
 */
public final class LongReadCanonicalAdapter {

    private static final String[] EVENT_MECHANISMS = {
            "TUTR",
            "ALE",
            "IPA",
            "internal_exon"
    };

    private LongReadCanonicalAdapter() {
    }

    public static void write(
            File outputDirectory,
            String sampleSheet,
            String gtf,
            String inputFormat,
            String control,
            String caseCondition,
            int clusterDistance,
            String clusterRepresentative,
            int minGeneReads,
            int minSiteReads,
            double minPau,
            double minDeltaPau,
            double fdr,
            double wpiThreshold,
            double internalExonThreshold,
            double aleThreshold,
            double ipaThreshold,
            double tutrThreshold,
            List<SiteCluster> retainedSites,
            DifferentialApa.Results differential,
            Map<String, ArchitectureAnalysis.Metrics> architecture
    ) throws IOException {
        Map<String, List<SiteCluster>> sitesByGene = groupSites(retainedSites);
        Map<String, List<DifferentialApa.SiteResult>> testsByGene =
                groupTests(differential.sites);

        List<Map<String, String>> apaRows =
                new ArrayList<Map<String, String>>();
        List<Map<String, String>> eventRows =
                new ArrayList<Map<String, String>>();
        List<Map<String, String>> siteRows =
                new ArrayList<Map<String, String>>();

        for (DifferentialApa.GeneResult gene : differential.genes) {
            List<SiteCluster> geneSites = sitesByGene.get(gene.geneId);
            if (geneSites == null || geneSites.isEmpty()) {
                continue;
            }

            List<DifferentialApa.SiteResult> geneTests =
                    testsByGene.get(gene.geneId);
            if (geneTests == null) {
                geneTests =
                        new ArrayList<DifferentialApa.SiteResult>();
            }

            ArchitectureAnalysis.Metrics metrics =
                    architecture.get(gene.geneId);

            apaRows.add(geneRow(
                    gene,
                    geneSites,
                    geneTests,
                    metrics
            ));

            for (String mechanism : EVENT_MECHANISMS) {
                List<SiteCluster> mechanismSites =
                        mechanismSites(geneSites, mechanism);

                if (mechanismSites.isEmpty()) {
                    continue;
                }

                eventRows.add(eventRow(
                        gene,
                        metrics,
                        mechanism,
                        mechanismSites,
                        geneTests
                ));
            }
        }

        for (DifferentialApa.SiteResult result : differential.sites) {
            siteRows.add(siteRow(result));
        }

        CanonicalOutputWriter.writeApa(
                new File(outputDirectory, "apa.tsv"),
                apaRows
        );
        CanonicalOutputWriter.writeEvents(
                new File(outputDirectory, "apa_events.tsv"),
                eventRows
        );
        CanonicalOutputWriter.writeSites(
                new File(outputDirectory, "apa_sites.tsv"),
                siteRows
        );

        LinkedHashMap<String, String> info =
                new LinkedHashMap<String, String>();

        info.put("apaura_version", Citation.VERSION);
        info.put(
                "mode",
                inputFormat != null
                        && inputFormat.toLowerCase(Locale.ROOT).startsWith("3seq:")
                        ? "run-3seq"
                        : ("isoquant".equalsIgnoreCase(inputFormat)
                                ? "run-isoquant"
                                : "run-direct")
        );
        info.put("input_format", inputFormat);
        info.put("sample_sheet", sampleSheet);
        info.put("annotation", gtf);
        info.put("control", control);
        info.put("case", caseCondition);
        info.put(
                "coordinate_system",
                "1-based genomic positions; start/end are analyzed PAS span, inclusive"
        );
        info.put("cluster_distance", Integer.toString(clusterDistance));
        info.put("cluster_representative", clusterRepresentative);
        info.put(
                "min_gene_reads_per_condition",
                Integer.toString(minGeneReads)
        );
        info.put(
                "min_site_reads_pooled",
                Integer.toString(minSiteReads)
        );
        info.put("min_pau_pooled", fmt(minPau));
        info.put("min_delta_pau", fmt(minDeltaPau));
        info.put("fdr_threshold", fmt(fdr));
        info.put("wpi_threshold", fmt(wpiThreshold));
        info.put(
                "min_internal_exon_delta",
                fmt(internalExonThreshold)
        );
        info.put("min_ale_switch", fmt(aleThreshold));
        info.put("min_ipa_delta", fmt(ipaThreshold));
        info.put("min_tutr_shift", fmt(tutrThreshold));
        info.put(
                "apa_summary_primary_mechanism",
                "among changed mechanisms: smallest mechanism FDR, "
                        + "then smallest mechanism p-value, "
                        + "then existing deterministic mechanism order"
        );
        info.put(
                "apa_summary_number_sites",
                "PAS assigned to primary_mechanism; NA when primary_mechanism=none"
        );
        info.put(
                "apa_summary_effect_size",
                "TUTR=tutrShiftMagnitude; ALE=aleSwitchIndex; "
                        + "IPA=abs(deltaIpaUsage); "
                        + "internal_exon=abs(deltaInternalUsage); "
                        + "NA when primary_mechanism=none"
        );
        info.put(
                "apa_summary_pvalue",
                inputFormat != null
                        && inputFormat.toLowerCase(Locale.ROOT).startsWith("3seq:")
                        ? "gene-level moderated Dirichlet-multinomial likelihood-ratio p-value"
                        : "minimum tested PAS p-value within primary_mechanism; "
                                + "NA when primary_mechanism=none"
        );
        info.put(
                "apa_summary_fdr",
                inputFormat != null
                        && inputFormat.toLowerCase(Locale.ROOT).startsWith("3seq:")
                        ? "Benjamini-Hochberg FDR across gene-level moderated Dirichlet-multinomial tests"
                        : "minimum BH-adjusted PAS FDR within primary_mechanism; "
                                + "NA when primary_mechanism=none"
        );
        info.put(
                "details_gene_effect_size",
                "overall gene APD_total_variation_distance"
        );
        info.put(
                "details_gene_pvalue",
                "tie-adjusted discrete two-sample ECDF maximum-distance test"
        );
        info.put(
                "details_gene_fdr",
                "Benjamini-Hochberg adjustment of gene-level ECDF test p-values"
        );
        info.put(
                "classification_support",
                "gene_ks_fdr<=threshold OR minimum_site_fdr<=threshold; "
                        + "mechanism calls additionally require existing mechanism effect thresholds"
        );
        info.put(
                "event_pvalue_semantics",
                "minimum tested site p-value within the mechanism"
        );
        info.put(
                "event_fdr_semantics",
                "minimum BH-adjusted site FDR within the mechanism"
        );
        info.put(
                "site_fdr",
                "Benjamini-Hochberg adjustment across tested PAS"
        );

        CanonicalOutputWriter.writeRunInfo(
                new File(outputDirectory, "run_info.txt"),
                info
        );
    }

    private static Map<String, String> geneRow(
            DifferentialApa.GeneResult gene,
            List<SiteCluster> sites,
            List<DifferentialApa.SiteResult> tests,
            ArchitectureAnalysis.Metrics metrics
    ) {
        Map<String, String> row =
                new LinkedHashMap<String, String>();

        String[] span = span(sites);
        SiteCluster first = sites.get(0);
        String primary = primaryMechanism(gene, tests);

        row.put("chrom", first.chrom);
        row.put("start", span[0]);
        row.put("end", span[1]);
        row.put("strand", Character.toString(first.strand.symbol()));
        row.put("gene_id", gene.geneId);
        row.put("gene_name", gene.geneName);
        row.put("biotype", gene.biotype);
        row.put("apa_architecture", gene.apaType);
        row.put("apa_change", gene.apaChange ? "yes" : "no");
        row.put("apa_direction", gene.apaDirection);
        row.put("changed_mechanisms", gene.changedMechanisms);
        row.put("primary_mechanism", primary);

        // RUN3SEQ_CANONICAL_DM
        if (Double.isFinite(gene.compositionPvalue)
                && Double.isFinite(gene.compositionFdr)) {
            row.put("number_sites", Integer.toString(gene.nSites));
            row.put("effect_size", fmt(gene.apd));
            row.put("pvalue", fmt(gene.compositionPvalue));
            row.put("fdr", fmt(gene.compositionFdr));
        } else if ("none".equals(primary)) {
            row.put("number_sites", "NA");
            row.put("effect_size", "NA");
            row.put("pvalue", "NA");
            row.put("fdr", "NA");
        } else {
            List<SiteCluster> primarySites =
                    mechanismSites(sites, primary);
            double primaryP =
                    minSiteStatistic(tests, primary, false);
            double primaryFdr =
                    minSiteStatistic(tests, primary, true);
            double classifiedFdr =
                    mechanismFdr(gene, primary);

            requireSameFiniteFdr(
                    gene.geneId,
                    primary,
                    classifiedFdr,
                    primaryFdr
            );

            row.put(
                    "number_sites",
                    Integer.toString(primarySites.size())
            );
            row.put(
                    "effect_size",
                    fmt(mechanismEffect(gene, metrics, primary))
            );
            row.put("pvalue", fmt(primaryP));
            row.put("fdr", fmt(primaryFdr));
        }

        row.put(
                "status",
                gene.nSites >= 2 ? "ok" : "insufficient_sites"
        );

        return row;
    }

    private static Map<String, String> eventRow(
            DifferentialApa.GeneResult gene,
            ArchitectureAnalysis.Metrics metrics,
            String mechanism,
            List<SiteCluster> sites,
            List<DifferentialApa.SiteResult> tests
    ) {
        Map<String, String> row =
                new LinkedHashMap<String, String>();

        String[] span = span(sites);
        SiteCluster first = sites.get(0);

        row.put("chrom", first.chrom);
        row.put("start", span[0]);
        row.put("end", span[1]);
        row.put("strand", Character.toString(first.strand.symbol()));
        row.put("gene_id", gene.geneId);
        row.put("gene_name", gene.geneName);
        row.put("biotype", gene.biotype);
        row.put("event_id", gene.geneId + "|" + mechanism);
        row.put("mechanism", mechanism);
        row.put("number_sites", Integer.toString(sites.size()));
        row.put(
                "apa_change",
                containsMechanism(gene.changedMechanisms, mechanism)
                        ? "yes"
                        : "no"
        );
        row.put("apa_direction", mechanismDirection(gene, mechanism));
        row.put("effect_size", fmt(mechanismEffect(gene, metrics, mechanism)));
        row.put(
                "pvalue",
                fmt(minSiteStatistic(tests, mechanism, false))
        );
        row.put(
                "fdr",
                fmt(minSiteStatistic(tests, mechanism, true))
        );
        row.put("status", mechanismStatus(tests, mechanism));

        return row;
    }

    private static Map<String, String> siteRow(
            DifferentialApa.SiteResult result
    ) {
        Map<String, String> row =
                new LinkedHashMap<String, String>();

        SiteCluster site = result.site;
        String mechanism = canonicalSiteMechanism(site);
        String eventId = canonicalSiteEventId(
                site.geneId,
                mechanism
        );

        row.put("chrom", site.chrom);
        row.put(
                "position",
                Integer.toString(site.representativePos0 + 1)
        );
        row.put("strand", Character.toString(site.strand.symbol()));
        row.put("gene_id", site.geneId);
        row.put("gene_name", site.geneName);
        row.put("event_id", eventId);
        row.put("mechanism", mechanism);
        row.put("feature_id", site.id);
        row.put("feature_rank", Integer.toString(site.siteRank));
        row.put("control_usage", fmt(result.controlPau));
        row.put("case_usage", fmt(result.casePau));
        row.put("delta_usage", fmt(result.delta));
        row.put("pvalue", fmt(result.p));
        row.put("fdr", fmt(result.q));
        row.put("status", result.status);

        return row;
    }

    private static Map<String, List<SiteCluster>> groupSites(
            List<SiteCluster> sites
    ) {
        Map<String, List<SiteCluster>> grouped =
                new LinkedHashMap<String, List<SiteCluster>>();

        for (SiteCluster site : sites) {
            if (".".equals(site.geneId)) {
                continue;
            }

            List<SiteCluster> values = grouped.get(site.geneId);
            if (values == null) {
                values = new ArrayList<SiteCluster>();
                grouped.put(site.geneId, values);
            }
            values.add(site);
        }

        return grouped;
    }

    private static Map<String, List<DifferentialApa.SiteResult>> groupTests(
            List<DifferentialApa.SiteResult> tests
    ) {
        Map<String, List<DifferentialApa.SiteResult>> grouped =
                new LinkedHashMap<String, List<DifferentialApa.SiteResult>>();

        for (DifferentialApa.SiteResult test : tests) {
            String geneId = test.site.geneId;
            List<DifferentialApa.SiteResult> values =
                    grouped.get(geneId);

            if (values == null) {
                values =
                        new ArrayList<DifferentialApa.SiteResult>();
                grouped.put(geneId, values);
            }

            values.add(test);
        }

        return grouped;
    }

    private static String canonicalSiteMechanism(
            SiteCluster site
    ) {
        String role = site.architectureRole;

        if (role == null || role.trim().isEmpty()) {
            return "other";
        }

        return role;
    }

    private static String canonicalSiteEventId(
            String geneId,
            String mechanism
    ) {
        for (String eventMechanism : EVENT_MECHANISMS) {
            if (eventMechanism.equals(mechanism)) {
                return geneId + "|" + mechanism;
            }
        }

        /*
         * ALE+TUTR contributes internally to both ALE and TUTR event
         * summaries, but represents one physical PAS in apa_sites.tsv.
         * Non-event roles likewise do not receive synthetic event IDs.
         */
        return "NA";
    }

    private static List<SiteCluster> mechanismSites(
            List<SiteCluster> sites,
            String mechanism
    ) {
        List<SiteCluster> result =
                new ArrayList<SiteCluster>();

        for (SiteCluster site : sites) {
            if (belongsToMechanism(site, mechanism)) {
                result.add(site);
            }
        }

        return result;
    }

    private static List<String> siteMechanisms(
            SiteCluster site
    ) {
        List<String> result = new ArrayList<String>();
        String role = site.architectureRole;

        if ("ALE+TUTR".equals(role)) {
            result.add("ALE");
            result.add("TUTR");
        } else if ("ALE".equals(role)) {
            result.add("ALE");
        } else if ("TUTR".equals(role)) {
            result.add("TUTR");
        } else if ("IPA".equals(role)) {
            result.add("IPA");
        } else if ("internal_exon".equals(role)) {
            result.add("internal_exon");
        } else if ("terminal".equals(role)) {
            result.add("terminal");
        } else if ("downstream".equals(role)) {
            result.add("downstream");
        } else {
            result.add("other");
        }

        return result;
    }

    private static boolean belongsToMechanism(
            SiteCluster site,
            String mechanism
    ) {
        String role = site.architectureRole;

        if ("TUTR".equals(mechanism)) {
            return "TUTR".equals(role)
                    || "ALE+TUTR".equals(role);
        }

        if ("ALE".equals(mechanism)) {
            return "ALE".equals(role)
                    || "ALE+TUTR".equals(role);
        }

        return mechanism.equals(role);
    }

    private static String primaryMechanism(
            DifferentialApa.GeneResult gene,
            List<DifferentialApa.SiteResult> tests
    ) {
        if (gene.changedMechanisms == null
                || "none".equals(gene.changedMechanisms)
                || gene.changedMechanisms.trim().isEmpty()) {
            return "none";
        }

        String[] mechanisms = gene.changedMechanisms.split("\\+");
        String best = mechanisms[0];
        double bestFdr = mechanismFdr(gene, best);
        double bestP = minSiteStatistic(
                tests,
                best,
                false
        );

        for (int index = 1; index < mechanisms.length; index++) {
            String candidate = mechanisms[index];
            double candidateFdr =
                    mechanismFdr(gene, candidate);
            double candidateP = minSiteStatistic(
                    tests,
                    candidate,
                    false
            );

            if (betterPrimary(
                    candidateFdr,
                    candidateP,
                    bestFdr,
                    bestP
            )) {
                best = candidate;
                bestFdr = candidateFdr;
                bestP = candidateP;
            }
        }

        return best;
    }

    private static boolean betterPrimary(
            double candidateFdr,
            double candidateP,
            double bestFdr,
            double bestP
    ) {
        if (finiteLess(candidateFdr, bestFdr)) {
            return true;
        }

        boolean fdrTie =
                (Double.isFinite(candidateFdr)
                        && Double.isFinite(bestFdr)
                        && Double.compare(candidateFdr, bestFdr) == 0)
                || (!Double.isFinite(candidateFdr)
                        && !Double.isFinite(bestFdr));

        return fdrTie && finiteLess(candidateP, bestP);
    }

    private static void requireSameFiniteFdr(
            String geneId,
            String mechanism,
            double classifiedFdr,
            double emittedFdr
    ) {
        if (!Double.isFinite(classifiedFdr)
                && !Double.isFinite(emittedFdr)) {
            return;
        }

        if (!Double.isFinite(classifiedFdr)
                || !Double.isFinite(emittedFdr)
                || Math.abs(classifiedFdr - emittedFdr) > 1.0e-12) {
            throw new IllegalStateException(
                    "Primary-mechanism FDR mismatch for "
                            + geneId
                            + " / "
                            + mechanism
                            + ": classifier="
                            + fmt(classifiedFdr)
                            + ", canonical="
                            + fmt(emittedFdr)
            );
        }
    }

    private static boolean finiteLess(
            double left,
            double right
    ) {
        if (Double.isFinite(left) && !Double.isFinite(right)) {
            return true;
        }

        return Double.isFinite(left)
                && Double.isFinite(right)
                && left < right;
    }

    private static boolean containsMechanism(
            String changedMechanisms,
            String mechanism
    ) {
        if (changedMechanisms == null
                || "none".equals(changedMechanisms)) {
            return false;
        }

        for (String value : changedMechanisms.split("\\+")) {
            if (mechanism.equals(value)) {
                return true;
            }
        }

        return false;
    }

    private static String mechanismDirection(
            DifferentialApa.GeneResult gene,
            String mechanism
    ) {
        if ("TUTR".equals(mechanism)) {
            return gene.tutrDirection;
        }
        if ("ALE".equals(mechanism)) {
            return gene.aleDirection;
        }
        if ("IPA".equals(mechanism)) {
            return gene.ipaDirection;
        }
        if ("internal_exon".equals(mechanism)) {
            return gene.internalExonDirection;
        }
        return "not_applicable";
    }

    private static double mechanismEffect(
            DifferentialApa.GeneResult gene,
            ArchitectureAnalysis.Metrics metrics,
            String mechanism
    ) {
        if ("TUTR".equals(mechanism)) {
            return gene.tutrShiftMagnitude;
        }

        if ("ALE".equals(mechanism)) {
            return gene.aleSwitchIndex;
        }

        if ("IPA".equals(mechanism)) {
            return Math.abs(gene.deltaIpaUsage);
        }

        if ("internal_exon".equals(mechanism)
                && metrics != null) {
            return Math.abs(metrics.deltaInternalUsage);
        }

        return Double.NaN;
    }

    private static double mechanismFdr(
            DifferentialApa.GeneResult gene,
            String mechanism
    ) {
        if ("TUTR".equals(mechanism)) {
            return gene.tutrMinFdr;
        }
        if ("ALE".equals(mechanism)) {
            return gene.aleMinFdr;
        }
        if ("IPA".equals(mechanism)) {
            return gene.ipaMinFdr;
        }
        if ("internal_exon".equals(mechanism)) {
            return gene.internalExonMinFdr;
        }
        return Double.NaN;
    }

    private static double minSiteStatistic(
            List<DifferentialApa.SiteResult> tests,
            String mechanism,
            boolean fdr
    ) {
        double minimum = Double.NaN;

        for (DifferentialApa.SiteResult test : tests) {
            if (!belongsToMechanism(test.site, mechanism)) {
                continue;
            }

            double value = fdr ? test.q : test.p;

            if (!Double.isFinite(value)) {
                continue;
            }

            if (!Double.isFinite(minimum) || value < minimum) {
                minimum = value;
            }
        }

        return minimum;
    }

    private static String mechanismStatus(
            List<DifferentialApa.SiteResult> tests,
            String mechanism
    ) {
        String fallback = "not_tested";

        for (DifferentialApa.SiteResult test : tests) {
            if (!belongsToMechanism(test.site, mechanism)) {
                continue;
            }

            if ("ok".equals(test.status)) {
                return "ok";
            }

            if ("insufficient_replicate_support".equals(test.status)) {
                fallback = "insufficient_replicate_support";
            } else if ("not_tested".equals(fallback)
                    && test.status != null
                    && !test.status.trim().isEmpty()) {
                fallback = test.status;
            }
        }

        return fallback;
    }

    private static String[] span(
            List<SiteCluster> sites
    ) {
        if (sites == null || sites.isEmpty()) {
            return new String[]{"NA", "NA"};
        }

        int minimum = Integer.MAX_VALUE;
        int maximum = Integer.MIN_VALUE;

        for (SiteCluster site : sites) {
            int position1 = site.representativePos0 + 1;
            minimum = Math.min(minimum, position1);
            maximum = Math.max(maximum, position1);
        }

        return new String[]{
                Integer.toString(minimum),
                Integer.toString(maximum)
        };
    }

    private static String fmt(double value) {
        if (!Double.isFinite(value)) {
            return "NA";
        }

        return String.format(
                Locale.ROOT,
                "%.10g",
                value
        );
    }
}
