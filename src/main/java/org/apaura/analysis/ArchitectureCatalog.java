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
        Map<String, List<TranscriptModel>> transcriptsByGene = groupByGene(transcripts);
        List<Entry> entries = new ArrayList<Entry>();

        for (List<TranscriptModel> geneTranscripts : transcriptsByGene.values()) {
            if (geneTranscripts.isEmpty()) {
                continue;
            }
            entries.add(buildEntry(geneTranscripts));
        }

        entries.sort((left, right) -> {
            int byChromosome = left.chrom.compareTo(right.chrom);
            return byChromosome != 0
                    ? byChromosome
                    : left.geneId.compareTo(right.geneId);
        });
        return entries;
    }

    private static Entry buildEntry(List<TranscriptModel> transcripts) {
        TranscriptModel first = transcripts.get(0);
        Entry entry = new Entry();
        entry.geneId = first.geneId;
        entry.geneName = first.geneName;
        entry.biotype = first.biotype;
        entry.chrom = first.chrom;
        entry.strand = first.strand;
        entry.transcripts = transcripts.size();

        Map<Integer, Set<Integer>> pasByTerminalBranch =
                new LinkedHashMap<Integer, Set<Integer>>();
        Set<Integer> allPas = new LinkedHashSet<Integer>();

        for (TranscriptModel transcript : transcripts) {
            pasByTerminalBranch.computeIfAbsent(
                    transcript.terminalExonAnchor0(),
                    key -> new LinkedHashSet<Integer>()
            ).add(transcript.pasPos0());
            allPas.add(transcript.pasPos0());
        }

        entry.terminalBranches = pasByTerminalBranch.size();
        entry.annotatedPas = allPas.size();
        for (Set<Integer> branchPas : pasByTerminalBranch.values()) {
            if (branchPas.size() >= 2) {
                entry.tutrBranches++;
            }
        }

        entry.ipaCompatibleEnds = countIpaCompatibleEnds(transcripts);
        setArchitecture(entry);
        return entry;
    }

    private static int countIpaCompatibleEnds(List<TranscriptModel> transcripts) {
        int count = 0;
        for (int index = 0; index < transcripts.size(); index++) {
            int pasPosition = transcripts.get(index).pasPos0();
            boolean compatible = false;

            for (int otherIndex = 0;
                    otherIndex < transcripts.size() && !compatible;
                    otherIndex++) {
                if (index == otherIndex) {
                    continue;
                }
                for (Interval intron : introns(transcripts.get(otherIndex))) {
                    if (intron.contains(pasPosition)) {
                        compatible = true;
                        break;
                    }
                }
            }
            if (compatible) {
                count++;
            }
        }
        return count;
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
