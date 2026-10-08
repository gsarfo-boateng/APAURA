package org.apaura.io;

import org.apaura.model.Interval;
import org.apaura.model.Strand;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public final class AlignmentRecord {
    public final String readName;
    public final String referenceName;
    public final int pos0;
    public final int mapq;
    public final int flag;
    public final int referenceConsumed;
    public final int queryLength;
    public final List<Interval> alignedBlocks;

    public AlignmentRecord(
            String readName,
            String referenceName,
            int pos0,
            int mapq,
            int flag,
            int referenceConsumed
    ) {
        this(
                readName,
                referenceName,
                pos0,
                mapq,
                flag,
                referenceConsumed,
                0,
                Collections.<Interval>emptyList()
        );
    }

    public AlignmentRecord(
            String readName,
            String referenceName,
            int pos0,
            int mapq,
            int flag,
            int referenceConsumed,
            int queryLength,
            List<Interval> alignedBlocks
    ) {
        this.readName = readName;
        this.referenceName = referenceName;
        this.pos0 = pos0;
        this.mapq = mapq;
        this.flag = flag;
        this.referenceConsumed = referenceConsumed;
        this.queryLength = queryLength;
        this.alignedBlocks = Collections.unmodifiableList(
                new ArrayList<Interval>(alignedBlocks)
        );
    }

    public boolean isUnmapped() {
        return (flag & 4) != 0 || referenceName == null || pos0 < 0;
    }

    public boolean isNegative() {
        return (flag & 16) != 0;
    }

    public boolean isPaired() {
        return (flag & 1) != 0;
    }

    public boolean isFirstOfPair() {
        return (flag & 64) != 0;
    }

    public boolean isSecondOfPair() {
        return (flag & 128) != 0;
    }

    public boolean isSecondary() {
        return (flag & 256) != 0;
    }

    public boolean isQcFail() {
        return (flag & 512) != 0;
    }

    public boolean isDuplicate() {
        return (flag & 1024) != 0;
    }

    public boolean isSupplementary() {
        return (flag & 2048) != 0;
    }

    public int end0Exclusive() {
        return pos0 + referenceConsumed;
    }

    public Strand alignmentStrand() {
        return isNegative() ? Strand.MINUS : Strand.PLUS;
    }
}
