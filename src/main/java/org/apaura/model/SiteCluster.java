package org.apaura.model;

import java.util.LinkedHashMap;
import java.util.Map;

public final class SiteCluster {
    public final String id;
    public final String chrom;
    public final Strand strand;
    public final int startPos0;
    public final int endPos0;
    public final int representativePos0;
    public final Map<String, Integer> counts = new LinkedHashMap<String, Integer>();

    public String geneId = ".";
    public String geneName = ".";
    public String biotype = ".";
    public String region = "INTERGENIC";
    public int siteRank = -1;
    public double normalizedPosition = Double.NaN;

    /** 5' anchor of the matched terminal exon, or -1 when not applicable. */
    /** Raw 5' anchor of the matched GTF terminal exon, or -1 when not applicable. */
    public int terminalExonAnchor0 = -1;

    /**
     * Reconciled terminal-structure anchor, or -1 when the site is not assigned
     * to a resolved terminal structure.
     */
    public int terminalStructureAnchor0 = -1;

    public boolean terminalExon = false;

    /** Gene-architecture role assigned after annotation. */
    public String architectureRole = "unclassified";

    /** 1-based terminal-exon branch rank in transcript orientation. */
    public int terminalBranchRank = -1;

    public SiteCluster(
            String id,
            String chrom,
            Strand strand,
            int startPos0,
            int endPos0,
            int representativePos0
    ) {
        this.id = id;
        this.chrom = chrom;
        this.strand = strand;
        this.startPos0 = startPos0;
        this.endPos0 = endPos0;
        this.representativePos0 = representativePos0;
    }

    public int totalCount() {
        int total = 0;
        for (int count : counts.values()) {
            total += count;
        }
        return total;
    }
}
