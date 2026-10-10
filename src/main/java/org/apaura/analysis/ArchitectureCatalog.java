package org.apaura.analysis;

import org.apaura.model.Interval;
import org.apaura.model.Strand;
import org.apaura.model.TranscriptModel;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/** Builds an annotation-derived inventory of potential APA architectures. */
public final class ArchitectureCatalog {

    private static final int DEFAULT_ANCHOR_TOLERANCE = 5;
    private static final int DEFAULT_PAS_MERGE_DISTANCE = 24;

    private ArchitectureCatalog() {
        // Utility class.
    }

    public static final class Entry {
        public String geneId;
        public String geneName;
        public String biotype;
        public String chrom;
        public String architectureClass;
        public String components;
        public Strand strand;
        public int geneStart1based;
        public int geneEnd1based;
        public int transcripts;
        public int terminalBranches;
        /** All unique annotated transcript-end PAS after PAS reconciliation. */
        public int annotatedPas;
        /** PAS retained in terminal TUTR/ALE architecture after IPA separation. */
        public int terminalPas;
        public int tutrBranches;
        public int ipaCompatibleEnds;
    }

    public static List<Entry> build(List<TranscriptModel> transcripts) {
        return build(
                transcripts,
                DEFAULT_ANCHOR_TOLERANCE,
                DEFAULT_PAS_MERGE_DISTANCE
        );
    }

    public static List<Entry> build(
            List<TranscriptModel> transcripts,
            int anchorTolerance,
            int pasMergeDistance
    ) {
        Map<String, List<TranscriptModel>> transcriptsByGene = groupByGene(transcripts);
        List<Entry> entries = new ArrayList<Entry>();

        for (List<TranscriptModel> geneTranscripts : transcriptsByGene.values()) {
            if (geneTranscripts.isEmpty()) {
                continue;
            }
            TerminalFamilyReconciler.Result allReconciliation =
                    TerminalFamilyReconciler.reconcile(
                            geneTranscripts,
                            anchorTolerance,
                            pasMergeDistance
                    );
            List<TranscriptModel> terminalTranscripts =
                    ShortReadIpa.terminalArchitectureTranscripts(geneTranscripts);
            TerminalFamilyReconciler.Result terminalReconciliation =
                    TerminalFamilyReconciler.reconcile(
                            terminalTranscripts,
                            anchorTolerance,
                            pasMergeDistance
                    );
            entries.add(buildEntry(
                    geneTranscripts,
                    terminalReconciliation.families,
                    allReconciliation.families
            ));
        }

        entries.sort((left, right) -> {
            int byChromosome = left.chrom.compareTo(right.chrom);
            return byChromosome != 0
                    ? byChromosome
                    : left.geneId.compareTo(right.geneId);
        });
        return entries;
    }

    private static Entry buildEntry(
            List<TranscriptModel> transcripts,
            List<TerminalFamilyReconciler.Family> terminalFamilies,
            List<TerminalFamilyReconciler.Family> allFamilies
    ) {
        TranscriptModel first = transcripts.get(0);
        Entry entry = new Entry();
        entry.geneId = first.geneId;
        entry.geneName = first.geneName;
        entry.biotype = preferredBiotype(transcripts);
        entry.chrom = first.chrom;
        entry.strand = first.strand;
        int geneStart0 = Integer.MAX_VALUE;
        int geneEnd0 = Integer.MIN_VALUE;
        for (TranscriptModel transcript : transcripts) {
            for (Interval exon : transcript.exons) {
                geneStart0 = Math.min(geneStart0, exon.start0);
                geneEnd0 = Math.max(geneEnd0, exon.end0);
            }
        }
        entry.geneStart1based = geneStart0 == Integer.MAX_VALUE ? 0 : geneStart0 + 1;
        entry.geneEnd1based = geneEnd0 == Integer.MIN_VALUE ? 0 : geneEnd0;
        entry.transcripts = transcripts.size();
        entry.terminalBranches = terminalFamilies.size();

        Set<String> terminalPas = physicalPas(terminalFamilies);
        for (TerminalFamilyReconciler.Family family : terminalFamilies) {
            if (family.pasPositions.size() >= 2) {
                entry.tutrBranches++;
            }
        }

        entry.annotatedPas = physicalPas(allFamilies).size();
        entry.terminalPas = terminalPas.size();
        entry.ipaCompatibleEnds = countIpaCompatibleEnds(transcripts);
        setArchitecture(entry);
        return entry;
    }

    private static Set<String> physicalPas(
            List<TerminalFamilyReconciler.Family> families
    ) {
        Set<String> positions = new LinkedHashSet<String>();
        for (TerminalFamilyReconciler.Family family : families) {
            for (Integer pas : family.pasPositions) {
                positions.add(
                        family.chrom + "\u0001" + family.strand.symbol() + "\u0001" + pas
                );
            }
        }
        return positions;
    }

    private static int countIpaCompatibleEnds(
            List<TranscriptModel> transcripts
    ) {
        Set<String> counted = new LinkedHashSet<String>();
        for (TranscriptModel transcript : transcripts) {
            if (ShortReadIpa.isStructuralIpaTranscript(transcript, transcripts)) {
                counted.add(
                        transcript.chrom
                                + "\u0001"
                                + transcript.strand.symbol()
                                + "\u0001"
                                + transcript.pasPos0()
                );
            }
        }
        return counted.size();
    }

    private static void setArchitecture(Entry entry) {
        boolean hasTutr = entry.tutrBranches > 0;
        boolean hasAle = entry.terminalBranches > 1;

        entry.components = ApaArchitecture.components(hasTutr, hasAle);
        entry.architectureClass = ApaArchitecture.classify(
                entry.terminalBranches,
                entry.tutrBranches,
                entry.terminalPas
        );
    }

    private static Map<String, List<TranscriptModel>> groupByGene(
            List<TranscriptModel> transcripts
    ) {
        Map<String, List<TranscriptModel>> grouped =
                new LinkedHashMap<String, List<TranscriptModel>>();
        for (TranscriptModel transcript : transcripts) {
            grouped.computeIfAbsent(
                    transcript.geneId,
                    key -> new ArrayList<TranscriptModel>()
            ).add(transcript);
        }
        return grouped;
    }

    private static String preferredBiotype(List<TranscriptModel> transcripts) {
        String biotype = ".";
        for (TranscriptModel transcript : transcripts) {
            if ("protein_coding".equals(transcript.biotype)) {
                return "protein_coding";
            }
            if (".".equals(biotype) && !".".equals(transcript.biotype)) {
                biotype = transcript.biotype;
            }
        }
        return biotype;
    }

    private static String join(List<String> values, String separator) {
        StringBuilder joined = new StringBuilder();
        for (String value : values) {
            if (joined.length() > 0) {
                joined.append(separator);
            }
            joined.append(value);
        }
        return joined.length() == 0 ? "none" : joined.toString();
    }
}