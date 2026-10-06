package org.apaura.model;

import java.util.ArrayList;
import java.util.List;

public final class GeneModel {
    public final String geneId;
    public String geneName;
    public String biotype;
    public final String chrom;
    public final Strand strand;
    public int start0 = Integer.MAX_VALUE;
    public int end0 = Integer.MIN_VALUE;
    public int cdsMin = Integer.MAX_VALUE;
    public int cdsMax = Integer.MIN_VALUE;
    public final List<Interval> exons = new ArrayList<Interval>();
    public final List<Interval> cds = new ArrayList<Interval>();
    public final List<Interval> utrs = new ArrayList<Interval>();
    public final List<Interval> explicitUtr3 = new ArrayList<Interval>();
    public final List<Interval> terminalExons = new ArrayList<Interval>();

    public GeneModel(
            String geneId,
            String geneName,
            String biotype,
            String chrom,
            Strand strand
    ) {
        this.geneId = geneId;
        this.geneName = geneName;
        this.biotype = biotype;
        this.chrom = chrom;
        this.strand = strand;
    }

    public void observe(int start, int end) {
        if (start < start0) {
            start0 = start;
        }
        if (end > end0) {
            end0 = end;
        }
    }

    public void addCds(int start, int end) {
        observe(start, end);
        cds.add(new Interval(start, end));
        cdsMin = Math.min(cdsMin, start);
        cdsMax = Math.max(cdsMax, end);
    }

    public boolean hasCds() {
        return cdsMin != Integer.MAX_VALUE;
    }

    public int transcriptThreePrimePos0() {
        return strand == Strand.MINUS ? start0 : end0 - 1;
    }

    public boolean overlaps(List<Interval> intervals, int position0) {
        for (Interval interval : intervals) {
            if (interval.contains(position0)) {
                return true;
            }
        }
        return false;
    }

    /** Returns true when a position lies in the 3' noncoding part of a terminal exon. */
    public boolean isThreePrimeUtrWithinTerminalExon(Interval terminalExon, int position0) {
        int localCdsStart = Integer.MAX_VALUE;
        int localCdsEnd = Integer.MIN_VALUE;
        boolean overlapsCds = false;

        for (Interval codingInterval : cds) {
            int overlapStart = Math.max(codingInterval.start0, terminalExon.start0);
            int overlapEnd = Math.min(codingInterval.end0, terminalExon.end0);
            if (overlapStart < overlapEnd) {
                overlapsCds = true;
                localCdsStart = Math.min(localCdsStart, overlapStart);
                localCdsEnd = Math.max(localCdsEnd, overlapEnd);
            }
        }

        if (!overlapsCds) {
            return true;
        }
        return strand == Strand.MINUS
                ? position0 < localCdsStart
                : position0 >= localCdsEnd;
    }

    public Interval bestTerminalExon(int position0) {
        Interval best = null;
        int bestDistance = Integer.MAX_VALUE;

        for (Interval terminalExon : terminalExons) {
            if (!terminalExon.contains(position0)) {
                continue;
            }

            int terminalThreePrimePosition = strand == Strand.MINUS
                    ? terminalExon.start0
                    : terminalExon.end0 - 1;
            int distance = Math.abs(position0 - terminalThreePrimePosition);
            if (distance < bestDistance) {
                best = terminalExon;
                bestDistance = distance;
            }
        }
        return best;
    }
}
