package org.apaura.model;

public final class EndKey implements Comparable<EndKey> {
    public final String chrom;
    public final int pos0;
    public final Strand strand;

    public EndKey(String chrom, int pos0, Strand strand) {
        this.chrom = chrom;
        this.pos0 = pos0;
        this.strand = strand;
    }

    @Override
    public int compareTo(EndKey other) {
        int chromosomeComparison = chrom.compareTo(other.chrom);
        if (chromosomeComparison != 0) {
            return chromosomeComparison;
        }

        int strandComparison = Character.compare(strand.symbol(), other.strand.symbol());
        if (strandComparison != 0) {
            return strandComparison;
        }

        return Integer.compare(pos0, other.pos0);
    }

    @Override
    public boolean equals(Object other) {
        if (this == other) {
            return true;
        }
        if (!(other instanceof EndKey)) {
            return false;
        }

        EndKey that = (EndKey) other;
        return pos0 == that.pos0
                && chrom.equals(that.chrom)
                && strand == that.strand;
    }

    @Override
    public int hashCode() {
        return 31 * (31 * chrom.hashCode() + pos0) + strand.hashCode();
    }
}
