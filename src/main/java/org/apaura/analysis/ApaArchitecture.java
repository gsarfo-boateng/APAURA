package org.apaura.analysis;

/**
 * Shared terminal APA architecture vocabulary and classification logic.
 *
 * Architecture describes structural organization of terminal poly(A) sites.
 * It is independent of case/control status, differential usage, p-values,
 * FDR, WPI, or other differential statistics.
 */
final class ApaArchitecture {

    static final String TUTR = "TUTR";
    static final String ALE = "ALE";
    static final String MIXED = "mixed";
    static final String NONE = "none";

    private ApaArchitecture() {
        // Utility class.
    }

    static String classify(boolean hasTutr, boolean hasAle) {
        if (hasTutr && hasAle) {
            return MIXED;
        }
        if (hasTutr) {
            return TUTR;
        }
        if (hasAle) {
            return ALE;
        }
        return NONE;
    }

    static String components(boolean hasTutr, boolean hasAle) {
        if (hasTutr && hasAle) {
            return "TUTR+ALE";
        }
        if (hasTutr) {
            return TUTR;
        }
        if (hasAle) {
            return ALE;
        }
        return NONE;
    }
}