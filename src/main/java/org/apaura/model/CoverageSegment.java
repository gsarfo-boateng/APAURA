package org.apaura.model;

public final class CoverageSegment {
    public final int index;
    public final int start0;
    public final int end0;
    public final String label;
    
    public CoverageSegment(int index, int start0, int end0, String label) {
        this.index = index;
        this.start0 = start0;
        this.end0 = end0;
        this.label = label;
    }
    
    public int length() {
        return end0 - start0;
    }
}