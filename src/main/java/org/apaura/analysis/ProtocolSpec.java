package org.apaura.analysis;

/** Endpoint interpretation used by APAURA run-3seq. */
public final class ProtocolSpec {

    public enum Boundary {
        TRANSCRIPT_3P,
        ALIGNMENT_5P
    }

    public final String name;
    public final Boundary boundary;
    public final boolean flipAlignmentStrand;
    public final boolean internalPrimingRecommended;
    public final int defaultClusterDistance;

    public ProtocolSpec(
            String name,
            Boundary boundary,
            boolean flipAlignmentStrand,
            boolean internalPrimingRecommended,
            int defaultClusterDistance
    ) {
        this.name = name;
        this.boundary = boundary;
        this.flipAlignmentStrand = flipAlignmentStrand;
        this.internalPrimingRecommended = internalPrimingRecommended;
        this.defaultClusterDistance = defaultClusterDistance;
    }
}
