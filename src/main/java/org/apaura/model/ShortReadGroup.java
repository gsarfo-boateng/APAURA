package org.apaura.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/** Tandem terminal-exon APA group with sites ordered from proximal to distal. */
public final class ShortReadGroup {
    public final String id;
    public final String geneId;
    public final String geneName;
    public final String biotype;
    public final String chrom;
    public final String scope;
    public final Strand strand;
    public final int anchor0;
    public final List<ShortReadSite> sites;
    public final List<CoverageSegment> segments;

    public ShortReadGroup(
            String id,
            String geneId,
            String geneName,
            String biotype,
            String chrom,
            Strand strand,
            String scope,
            int anchor0,
            List<ShortReadSite> sites,
            List<CoverageSegment> segments
    ) {
        this.id = id;
        this.geneId = geneId;
        this.geneName = geneName;
        this.biotype = biotype;
        this.chrom = chrom;
        this.strand = strand;
        this.scope = scope;
        this.anchor0 = anchor0;
        this.sites = Collections.unmodifiableList(new ArrayList<ShortReadSite>(sites));
        this.segments = Collections.unmodifiableList(new ArrayList<CoverageSegment>(segments));
    }
}
