package org.apaura.annotation;

import org.apaura.model.GeneModel;
import org.apaura.model.Interval;
import org.apaura.model.SiteCluster;
import org.apaura.model.Strand;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/** GTF-backed annotation of PAS coordinates to genes and transcript regions. */
public final class GtfAnnotation {
    private final Map<String, List<GeneModel>> genesByChromosome =
            new HashMap<String, List<GeneModel>>();

    private static final class TranscriptBuffer {
        final String geneId;
        final Strand strand;
        final List<Interval> exons = new ArrayList<Interval>();

        TranscriptBuffer(String geneId, Strand strand) {
            this.geneId = geneId;
            this.strand = strand;
        }

        Interval terminalExon() {
            if (exons.isEmpty()) {
                return null;
            }

            Interval terminal = exons.get(0);
            for (Interval exon : exons) {
                if (strand == Strand.MINUS) {
                    if (exon.start0 < terminal.start0) {
                        terminal = exon;
                    }
                } else if (exon.end0 > terminal.end0) {
                    terminal = exon;
                }
            }
            return terminal;
        }
    }

    private static final class Candidate {
        final GeneModel gene;
        final String region;
        final int priority;
        final int distance;
        final boolean terminalExon;
        final int terminalExonAnchor0;

        Candidate(
                GeneModel gene,
                String region,
                int priority,
                int distance,
                boolean terminalExon,
                int terminalExonAnchor0
        ) {
            this.gene = gene;
            this.region = region;
            this.priority = priority;
            this.distance = distance;
            this.terminalExon = terminalExon;
            this.terminalExonAnchor0 = terminalExonAnchor0;
        }
    }

    private GtfAnnotation() {
    }

    public static GtfAnnotation load(String path) throws IOException {
        GtfAnnotation annotation = new GtfAnnotation();
        Map<String, GeneModel> genes = new LinkedHashMap<String, GeneModel>();
        Map<String, TranscriptBuffer> transcripts =
                new LinkedHashMap<String, TranscriptBuffer>();

        try (BufferedReader reader = openReader(path)) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '#') {
                    continue;
                }
                parseRecord(line, genes, transcripts);
            }
        }

        attachTranscriptTerminalExons(genes, transcripts);
        addFallbackTerminalExons(genes);
        annotation.indexGenes(genes);

        System.err.println(
                "[APAURA] loaded GTF genes="
                        + genes.size()
                        + " transcript_terminal_exons="
                        + transcripts.size()
        );
        return annotation;
    }

    private static BufferedReader openReader(String path) throws IOException {
        InputStream input = new FileInputStream(path);
        if (path.endsWith(".gz")) {
            input = new GZIPInputStream(input);
        }
        return new BufferedReader(new InputStreamReader(input, StandardCharsets.UTF_8));
    }

    private static void parseRecord(
            String line,
            Map<String, GeneModel> genes,
            Map<String, TranscriptBuffer> transcripts
    ) {
        String[] fields = line.split("\t");
        if (fields.length < 9) {
            return;
        }

        String feature = fields[2].toLowerCase(Locale.ROOT);
        int start0 = Integer.parseInt(fields[3]) - 1;
        int end0 = Integer.parseInt(fields[4]);
        Strand strand = Strand.parse(fields[6]);
        Map<String, String> attributes = parseAttributes(fields[8]);
        String geneId = attributes.get("gene_id");
        if (geneId == null) {
            return;
        }

        String geneName = attributes.containsKey("gene_name")
                ? attributes.get("gene_name")
                : geneId;
        String biotype = firstPresent(
                attributes,
                "gene_type",
                "gene_biotype",
                "transcript_type",
                "transcript_biotype"
        );

        GeneModel gene = genes.get(geneId);
        if (gene == null) {
            gene = new GeneModel(
                    geneId,
                    geneName,
                    biotype == null ? "." : biotype,
                    fields[0],
                    strand
            );
            genes.put(geneId, gene);
        }

        gene.observe(start0, end0);
        updateGeneMetadata(gene, geneName, biotype);
        addFeature(feature, start0, end0, strand, attributes, gene, transcripts);
    }

    private static void updateGeneMetadata(
            GeneModel gene,
            String geneName,
            String biotype
    ) {
        if (gene.geneName == null || gene.geneName.equals(gene.geneId)) {
            gene.geneName = geneName;
        }
        if ((gene.biotype == null || gene.biotype.equals(".")) && biotype != null) {
            gene.biotype = biotype;
        }
    }

    private static void addFeature(
            String feature,
            int start0,
            int end0,
            Strand strand,
            Map<String, String> attributes,
            GeneModel gene,
            Map<String, TranscriptBuffer> transcripts
    ) {
        if (feature.equals("exon")) {
            Interval exon = new Interval(start0, end0);
            gene.exons.add(exon);

            String transcriptId = attributes.get("transcript_id");
            if (transcriptId != null) {
                String transcriptKey = gene.geneId + "\u0001" + transcriptId;
                TranscriptBuffer transcript = transcripts.get(transcriptKey);
                if (transcript == null) {
                    transcript = new TranscriptBuffer(gene.geneId, strand);
                    transcripts.put(transcriptKey, transcript);
                }
                transcript.exons.add(exon);
            }
            return;
        }

        if (feature.equals("cds")) {
            gene.addCds(start0, end0);
        } else if (feature.equals("utr")) {
            gene.utrs.add(new Interval(start0, end0));
        } else if (feature.equals("three_prime_utr")
                || feature.equals("3utr")
                || feature.equals("3'utr")) {
            gene.explicitUtr3.add(new Interval(start0, end0));
        }
    }

    private static void attachTranscriptTerminalExons(
            Map<String, GeneModel> genes,
            Map<String, TranscriptBuffer> transcripts
    ) {
        for (TranscriptBuffer transcript : transcripts.values()) {
            GeneModel gene = genes.get(transcript.geneId);
            Interval terminalExon = transcript.terminalExon();
            if (gene != null && terminalExon != null) {
                addUnique(gene.terminalExons, terminalExon);
            }
        }
    }

    /** Fallback for GTF files that do not provide transcript identifiers. */
    private static void addFallbackTerminalExons(Map<String, GeneModel> genes) {
        for (GeneModel gene : genes.values()) {
            if (!gene.terminalExons.isEmpty() || gene.exons.isEmpty()) {
                continue;
            }

            int extreme = gene.strand == Strand.MINUS
                    ? Integer.MAX_VALUE
                    : Integer.MIN_VALUE;

            for (Interval exon : gene.exons) {
                extreme = gene.strand == Strand.MINUS
                        ? Math.min(extreme, exon.start0)
                        : Math.max(extreme, exon.end0);
            }

            for (Interval exon : gene.exons) {
                boolean isTerminal = gene.strand == Strand.MINUS
                        ? exon.start0 == extreme
                        : exon.end0 == extreme;
                if (isTerminal) {
                    addUnique(gene.terminalExons, exon);
                }
            }
        }
    }

    private void indexGenes(Map<String, GeneModel> genes) {
        for (GeneModel gene : genes.values()) {
            List<GeneModel> chromosomeGenes = genesByChromosome.get(gene.chrom);
            if (chromosomeGenes == null) {
                chromosomeGenes = new ArrayList<GeneModel>();
                genesByChromosome.put(gene.chrom, chromosomeGenes);
            }
            chromosomeGenes.add(gene);
        }

        for (List<GeneModel> chromosomeGenes : genesByChromosome.values()) {
            Collections.sort(chromosomeGenes, new Comparator<GeneModel>() {
                @Override
                public int compare(GeneModel left, GeneModel right) {
                    return Integer.compare(left.start0, right.start0);
                }
            });
        }
    }

    private static void addUnique(List<Interval> intervals, Interval candidate) {
        for (Interval interval : intervals) {
            if (interval.start0 == candidate.start0 && interval.end0 == candidate.end0) {
                return;
            }
        }
        intervals.add(candidate);
    }

    private static String firstPresent(
            Map<String, String> attributes,
            String... keys
    ) {
        for (String key : keys) {
            if (attributes.containsKey(key)) {
                return attributes.get(key);
            }
        }
        return null;
    }

    private static Map<String, String> parseAttributes(String attributesText) {
        Map<String, String> attributes = new HashMap<String, String>();

        for (String rawAttribute : attributesText.split(";")) {
            String attribute = rawAttribute.trim();
            if (attribute.isEmpty()) {
                continue;
            }

            int space = attribute.indexOf(' ');
            int equals = attribute.indexOf('=');
            int separator = (space >= 0 && (equals < 0 || space < equals))
                    ? space
                    : equals;
            if (separator < 0) {
                continue;
            }

            String key = attribute.substring(0, separator).trim();
            String value = attribute.substring(separator + 1).trim();
            if (value.startsWith("\"")
                    && value.endsWith("\"")
                    && value.length() >= 2) {
                value = value.substring(1, value.length() - 1);
            }
            attributes.put(key, value);
        }
        return attributes;
    }

    public void annotate(List<SiteCluster> sites, int slop) {
        for (SiteCluster site : sites) {
            List<GeneModel> chromosomeGenes = genesByChromosome.get(site.chrom);
            if (chromosomeGenes == null) {
                continue;
            }

            Candidate best = bestCandidate(site, chromosomeGenes, slop);
            if (best != null) {
                applyCandidate(site, best);
            }
        }
        assignRanks(sites);
    }

    private static Candidate bestCandidate(
            SiteCluster site,
            List<GeneModel> genes,
            int slop
    ) {
        Candidate best = null;
        int position0 = site.representativePos0;
        int upperBound = upperBoundStart(genes, position0 + slop);

        for (int i = upperBound - 1; i >= 0; i--) {
            GeneModel gene = genes.get(i);
            if (gene.end0 + slop < position0) {
                break;
            }
            if (!strandCompatible(site.strand, gene.strand)) {
                continue;
            }

            Candidate candidate = classify(gene, position0, slop);
            if (isBetterCandidate(candidate, best)) {
                best = candidate;
            }
        }
        return best;
    }

    private static boolean strandCompatible(Strand siteStrand, Strand geneStrand) {
        return geneStrand == Strand.UNKNOWN
                || siteStrand == Strand.UNKNOWN
                || geneStrand == siteStrand;
    }

    private static boolean isBetterCandidate(Candidate candidate, Candidate currentBest) {
        if (candidate == null) {
            return false;
        }
        return currentBest == null
                || candidate.priority > currentBest.priority
                || (candidate.priority == currentBest.priority
                && candidate.distance < currentBest.distance);
    }

    private static void applyCandidate(SiteCluster site, Candidate candidate) {
        site.geneId = candidate.gene.geneId;
        site.geneName = candidate.gene.geneName;
        site.biotype = candidate.gene.biotype;
        site.region = candidate.region;
        site.terminalExon = candidate.terminalExon;
        site.terminalExonAnchor0 = candidate.terminalExonAnchor0;
    }

    private static int upperBoundStart(List<GeneModel> genes, int position0) {
        int lower = 0;
        int upper = genes.size();

        while (lower < upper) {
            int middle = (lower + upper) >>> 1;
            if (genes.get(middle).start0 <= position0) {
                lower = middle + 1;
            } else {
                upper = middle;
            }
        }
        return lower;
    }

    private static Candidate classify(GeneModel gene, int position0, int slop) {
        int distance = Math.abs(position0 - gene.transcriptThreePrimePos0());
        Interval terminalExon = gene.bestTerminalExon(position0);
        boolean isTerminalExon = terminalExon != null;
        int terminalAnchor0 = isTerminalExon
                ? terminalAnchor(gene, terminalExon)
                : -1;

        if (gene.overlaps(gene.explicitUtr3, position0)) {
            return candidate(
                    gene,
                    "THREE_PRIME_UTR",
                    60,
                    distance,
                    isTerminalExon,
                    terminalAnchor0
            );
        }

        if (gene.overlaps(gene.utrs, position0)
                && gene.hasCds()
                && isThreePrimeUtrByCdsBoundary(gene, position0)) {
            return candidate(
                    gene,
                    "THREE_PRIME_UTR",
                    58,
                    distance,
                    isTerminalExon,
                    terminalAnchor0
            );
        }

        if (isTerminalExon
                && gene.hasCds()
                && gene.isThreePrimeUtrWithinTerminalExon(terminalExon, position0)) {
            return candidate(
                    gene,
                    "THREE_PRIME_UTR",
                    56,
                    distance,
                    true,
                    terminalAnchor0
            );
        }

        if (isTerminalExon && !gene.hasCds()) {
            return candidate(
                    gene,
                    "TERMINAL_EXON_NONCODING",
                    52,
                    distance,
                    true,
                    terminalAnchor0
            );
        }

        if (gene.overlaps(gene.exons, position0)) {
            return candidate(
                    gene,
                    gene.hasCds() ? "EXON_CODING" : "EXON_NONCODING",
                    40,
                    distance,
                    false,
                    -1
            );
        }

        if (position0 >= gene.start0 && position0 < gene.end0) {
            return candidate(gene, "INTRON", 30, distance, false, -1);
        }

        if (isDownstreamThreePrime(gene, position0, slop)) {
            return candidate(gene, "DOWNSTREAM_3P", 20, distance, false, -1);
        }

        return null;
    }

    private static Candidate candidate(
            GeneModel gene,
            String region,
            int priority,
            int distance,
            boolean terminalExon,
            int terminalAnchor0
    ) {
        return new Candidate(
                gene,
                region,
                priority,
                distance,
                terminalExon,
                terminalAnchor0
        );
    }

    private static int terminalAnchor(GeneModel gene, Interval terminalExon) {
        return gene.strand == Strand.MINUS
                ? terminalExon.end0
                : terminalExon.start0;
    }

    private static boolean isThreePrimeUtrByCdsBoundary(
            GeneModel gene,
            int position0
    ) {
        return gene.strand == Strand.MINUS
                ? position0 < gene.cdsMin
                : position0 >= gene.cdsMax;
    }

    private static boolean isDownstreamThreePrime(
            GeneModel gene,
            int position0,
            int slop
    ) {
        if (gene.strand == Strand.MINUS) {
            return position0 < gene.start0 && gene.start0 - position0 <= slop;
        }
        return position0 >= gene.end0 && position0 - gene.end0 < slop;
    }

    private static void assignRanks(List<SiteCluster> sites) {
        Map<String, List<SiteCluster>> sitesByGene =
                new HashMap<String, List<SiteCluster>>();

        for (SiteCluster site : sites) {
            if (".".equals(site.geneId)) {
                continue;
            }
            List<SiteCluster> geneSites = sitesByGene.get(site.geneId);
            if (geneSites == null) {
                geneSites = new ArrayList<SiteCluster>();
                sitesByGene.put(site.geneId, geneSites);
            }
            geneSites.add(site);
        }

        for (List<SiteCluster> geneSites : sitesByGene.values()) {
            sortInTranscriptOrder(geneSites);
            int siteCount = geneSites.size();

            for (int i = 0; i < siteCount; i++) {
                SiteCluster site = geneSites.get(i);
                site.siteRank = i + 1;
                site.normalizedPosition = siteCount == 1
                        ? 0.5
                        : (double) i / (siteCount - 1);
            }
        }
    }

    private static void sortInTranscriptOrder(List<SiteCluster> sites) {
        Collections.sort(sites, new Comparator<SiteCluster>() {
            @Override
            public int compare(SiteCluster left, SiteCluster right) {
                int comparison = Integer.compare(
                        left.representativePos0,
                        right.representativePos0
                );
                return left.strand == Strand.MINUS ? -comparison : comparison;
            }
        });
    }
}
