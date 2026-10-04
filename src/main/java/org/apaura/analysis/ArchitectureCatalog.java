package org.apaura.analysis;

import org.apaura.model.Interval;
import org.apaura.model.Strand;
import org.apaura.model.TranscriptModel;

import java.util.ArrayList;
import java.util.Comparator;
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
        public int transcripts;
        public int terminalBranches;
        public int annotatedPas;
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
            TerminalFamilyReconciler.Result reconciliation =
                    TerminalFamilyReconciler.reconcile(
                            geneTranscripts,
                            anchorTolerance,
                            pasMergeDistance
                    );
            entries.add(buildEntry(geneTranscripts, reconciliation.families));
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
            List<TerminalFamilyReconciler.Family> families
    ) {
        TranscriptModel first = transcripts.get(0);
        Entry entry = new Entry();
        entry.geneId = first.geneId;
        entry.geneName = first.geneName;
        entry.biotype = preferredBiotype(transcripts);
        entry.chrom = first.chrom;
        entry.strand = first.strand;
        entry.transcripts = transcripts.size();
        entry.terminalBranches = families.size();

        Set<String> physicalPas = new LinkedHashSet<String>();
        for (TerminalFamilyReconciler.Family family : families) {
            if (family.pasPositions.size() >= 2) {
                entry.tutrBranches++;
            }
            for (Integer pas : family.pasPositions) {
                physicalPas.add(
                        family.chrom + "\u0001" + family.strand.symbol() + "\u0001" + pas
                );
            }
        }

        entry.annotatedPas = physicalPas.size();
        entry.ipaCompatibleEnds = countIpaCompatibleEnds(families, transcripts);
        setArchitecture(entry);
        return entry;
    }

    private static int countIpaCompatibleEnds(
            List<TerminalFamilyReconciler.Family> families,
            List<TranscriptModel> transcripts
    ) {
        Set<String> counted = new LinkedHashSet<String>();
        for (TerminalFamilyReconciler.Family family : families) {
            for (Integer pas : family.pasPositions) {
                if (fallsInAnyIntron(pas, family.chrom, family.strand, transcripts)) {
                    counted.add(
                            family.chrom
                                    + "\u0001"
                                    + family.strand.symbol()
                                    + "\u0001"
                                    + pas
                    );
                }
            }
        }
        return counted.size();
    }

    private static boolean fallsInAnyIntron(
            int pasPosition,
            String chrom,
            Strand strand,
            List<TranscriptModel> transcripts
    ) {
        for (TranscriptModel transcript : transcripts) {
            if (!chrom.equals(transcript.chrom) || strand != transcript.strand) {
                continue;
            }
            for (Interval intron : introns(transcript)) {
                if (intron.contains(pasPosition)) {
                    return true;
                }
            }
        }
        return false;
    }

    private static void setArchitecture(Entry entry) {
        boolean hasTutr = entry.tutrBranches > 0;
        boolean hasAle = entry.terminalBranches > 1;

        List<String> components = new ArrayList<String>();
        if (hasTutr) {
            components.add("TUTR");
        }
        if (hasAle) {
            components.add("ALE");
        }

        entry.components = join(components, "+");
        if (hasTutr && hasAle) {
            entry.architectureClass = "mixed";
        } else if (hasTutr) {
            entry.architectureClass = "TUTR";
        } else if (hasAle) {
            entry.architectureClass = "ALE";
        } else {
            entry.architectureClass = "none";
        }
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

    private static List<Interval> introns(TranscriptModel transcript) {
        List<Interval> exons = new ArrayList<Interval>(transcript.exons);
        exons.sort(Comparator.comparingInt(interval -> interval.start0));

        List<Interval> introns = new ArrayList<Interval>();
        for (int index = 1; index < exons.size(); index++) {
            Interval previous = exons.get(index - 1);
            Interval current = exons.get(index);
            if (previous.end0 < current.start0) {
                introns.add(new Interval(previous.end0, current.start0));
            }
        }
        return introns;
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