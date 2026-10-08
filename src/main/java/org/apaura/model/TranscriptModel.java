package org.apaura.model;

import java.util.ArrayList;
import java.util.List;

/** Transcript structure parsed from a GTF for annotated APA analysis. */
public final class TranscriptModel {
    public final String transcriptId;
    public final String geneId;
    public final String geneName;
    public final String biotype;
    public final String chrom;
    public final Strand strand;
    public final List<Interval> exons = new ArrayList<Interval>();
    public int cdsMin = Integer.MAX_VALUE;
    public int cdsMax = Integer.MIN_VALUE;

    public TranscriptModel(
            String transcriptId,
            String geneId,
            String geneName,
            String biotype,
            String chrom,
            Strand strand
    ) {
        this.transcriptId = transcriptId;
        this.geneId = geneId;
        this.geneName = geneName;
        this.biotype = biotype;
        this.chrom = chrom;
        this.strand = strand;
    }

    public void addExon(int start0, int end0) {
        exons.add(new Interval(start0, end0));
    }

    public void addCds(int start0, int end0) {
        cdsMin = Math.min(cdsMin, start0);
        cdsMax = Math.max(cdsMax, end0);
    }

    public boolean hasCds() {
        return cdsMin != Integer.MAX_VALUE;
    }

    public Interval terminalExon() {
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

    /** Single-base transcript end coordinate in 0-based genomic coordinates. */
    public int pasPos0() {
        Interval terminal = terminalExon();
        return strand == Strand.MINUS ? terminal.start0 : terminal.end0 - 1;
    }

    /** 5' boundary of the terminal exon as a half-open genomic boundary. */
    public int terminalExonAnchor0() {
        Interval terminal = terminalExon();
        return strand == Strand.MINUS ? terminal.end0 : terminal.start0;
    }
}
