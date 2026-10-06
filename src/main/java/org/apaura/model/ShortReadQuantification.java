package org.apaura.model;

public final class ShortReadQuantification {
    public final ShortReadGroup group;
    public final long[] alignedBases;
    public final double[] rawCoverage;
    public final double[] fittedCoverage;
    public final double[] siteAbundance;
    public final double[] sitePau;
    public final int[] pseudoCounts;
    public final int commonEvidenceReads;

    public ShortReadQuantification(
            ShortReadGroup group,
            long[] alignedBases,
            double[] rawCoverage,
            double[] fittedCoverage,
            double[] siteAbundance,
            double[] sitePau,
            int[] pseudoCounts,
            int commonEvidenceReads
    ) {
        this.group = group;
        this.alignedBases = alignedBases;
        this.rawCoverage = rawCoverage;
        this.fittedCoverage = fittedCoverage;
        this.siteAbundance = siteAbundance;
        this.sitePau = sitePau;
        this.pseudoCounts = pseudoCounts;
        this.commonEvidenceReads = commonEvidenceReads;
    }
}
