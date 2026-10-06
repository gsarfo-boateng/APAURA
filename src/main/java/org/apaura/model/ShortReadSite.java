package org.apaura.model;

public final class ShortReadSite {
    public final String id;
    public final int pos0;
    public final int rank;
    public final double normalizedPosition;

    public ShortReadSite(String id, int pos0, int rank, double normalizedPosition) {
        this.id = id;
        this.pos0 = pos0;
        this.rank = rank;
        this.normalizedPosition = normalizedPosition;
    }
}