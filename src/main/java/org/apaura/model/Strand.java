package org.apaura.model;

public enum Strand {
    PLUS('+'),
    MINUS('-'),
    UNKNOWN('.');

    private final char symbol;

    Strand(char symbol) {
        this.symbol = symbol;
    }

    public char symbol() {
        return symbol;
    }

    public static Strand fromFlag(boolean negative) {
        return negative ? MINUS : PLUS;
    }

    public static Strand parse(String value) {
        if ("+".equals(value) || "plus".equalsIgnoreCase(value)) {
            return PLUS;
        }
        if ("-".equals(value) || "minus".equalsIgnoreCase(value)) {
            return MINUS;
        }
        return UNKNOWN;
    }

    public Strand flip() {
        if (this == PLUS) {
            return MINUS;
        }
        if (this == MINUS) {
            return PLUS;
        }
        return UNKNOWN;
    }
}
