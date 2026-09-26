package org.apaura.model;

public final class Interval {
    public final int start0;
    public final int end0;
    
    public Interval(int start0, int end0) {
        this.start0 = start0;
        this.end0 = end0;
    }
    
    public boolean contains(int p) {
        return p >= start0 && p < end0;
    }
}