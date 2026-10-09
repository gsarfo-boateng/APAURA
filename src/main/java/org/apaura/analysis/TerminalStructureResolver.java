package org.apaura.analysis;

import org.apaura.model.SiteCluster;
import org.apaura.model.TranscriptModel;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Maps raw GTF terminal-exon anchors on observed PAS to reconciled
 * terminal structures.
 */
final class TerminalStructureResolver {

    static final int DEFAULT_ANCHOR_TOLERANCE = 5;
    static final int DEFAULT_PAS_MERGE_DISTANCE = 24;

    private TerminalStructureResolver() {
        // Utility class.
    }

    static void assign(
            List<TranscriptModel> transcripts,
            List<SiteCluster> sites
    ) {
        assign(
                transcripts,
                sites,
                DEFAULT_ANCHOR_TOLERANCE,
                DEFAULT_PAS_MERGE_DISTANCE
        );
    }

    static void assign(
            List<TranscriptModel> transcripts,
            List<SiteCluster> sites,
            int anchorTolerance,
            int pasMergeDistance
    ) {
        List<TranscriptModel> terminalTranscripts =
                ShortReadIpa.terminalArchitectureTranscripts(
                        transcripts
                );

        Map<String, List<TranscriptModel>> transcriptsByGene =
                groupTranscriptsByGene(terminalTranscripts);

        Map<String, Map<Integer, Integer>> structureByRawAnchor =
                new LinkedHashMap<String, Map<Integer, Integer>>();

        for (Map.Entry<String, List<TranscriptModel>> entry
                : transcriptsByGene.entrySet()) {

            TerminalFamilyReconciler.Result reconciliation =
                    TerminalFamilyReconciler.reconcile(
                            entry.getValue(),
                            anchorTolerance,
                            pasMergeDistance
                    );

            Map<Integer, Integer> geneStructures =
                    new LinkedHashMap<Integer, Integer>();

            for (TerminalFamilyReconciler.Family family
                    : reconciliation.families) {

                for (Integer rawAnchor : family.sourceAnchors) {
                    geneStructures.put(
                            rawAnchor,
                            family.anchor0
                    );
                }
            }

            structureByRawAnchor.put(
                    entry.getKey(),
                    geneStructures
            );
        }

        for (SiteCluster site : sites) {
            site.terminalStructureAnchor0 = -1;

            if (!site.terminalExon
                    || site.terminalExonAnchor0 < 0
                    || ".".equals(site.geneId)) {
                continue;
            }

            Map<Integer, Integer> geneStructures =
                    structureByRawAnchor.get(site.geneId);

            if (geneStructures == null) {
                continue;
            }

            Integer structureAnchor =
                    geneStructures.get(
                            site.terminalExonAnchor0
                    );

            if (structureAnchor != null) {
                site.terminalStructureAnchor0 =
                        structureAnchor;
            }
        }
    }

    private static Map<String, List<TranscriptModel>>
    groupTranscriptsByGene(
            List<TranscriptModel> transcripts
    ) {
        Map<String, List<TranscriptModel>> grouped =
                new LinkedHashMap<String, List<TranscriptModel>>();

        for (TranscriptModel transcript : transcripts) {
            List<TranscriptModel> geneTranscripts =
                    grouped.get(transcript.geneId);

            if (geneTranscripts == null) {
                geneTranscripts =
                        new ArrayList<TranscriptModel>();

                grouped.put(
                        transcript.geneId,
                        geneTranscripts
                );
            }

            geneTranscripts.add(transcript);
        }

        return grouped;
    }
}