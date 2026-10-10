package org.apaura.analysis;

/**
 * Shared terminal APA architecture vocabulary and classification logic.
 *
 * Architecture describes structural organization of terminal poly(A) sites.
 * It is independent of case/control status, differential usage, p-values,
 * FDR, WPI, or other differential statistics.
 */
final class ApaArchitecture {

    static final String SINGLE_PAS = "SINGLE_PAS";
    static final String TUTR = "TUTR";
    static final String ALE = "ALE";
    static final String MIXED = "mixed";
    static final String UNRESOLVED = "UNRESOLVED";

    private ApaArchitecture() {
        // Utility class.
    }

    /**
     * Classifies terminal APA structure from resolved terminal structures
     * and their retained terminal PAS.
     *
     * @param terminalStructures number of resolved terminal structures
     * @param tutrStructures number of terminal structures containing >=2 PAS
     * @param terminalPas number of retained terminal PAS
     */
    static String classify(
            int terminalStructures,
            int tutrStructures,
            int terminalPas
    ) {
        if (terminalStructures <= 0 || terminalPas <= 0) {
            return UNRESOLVED;
        }

        boolean hasTutr = tutrStructures > 0;
        boolean hasAle = terminalStructures > 1;

        if (hasTutr && hasAle) {
            return MIXED;
        }

        if (hasTutr) {
            return TUTR;
        }

        if (hasAle) {
            return ALE;
        }

        if (terminalStructures == 1 && terminalPas == 1) {
            return SINGLE_PAS;
        }

        return UNRESOLVED;
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
        return "none";
    }

    static boolean isAlternativeTerminalArchitecture(String architectureClass) {
        return TUTR.equals(architectureClass)
                || ALE.equals(architectureClass)
                || MIXED.equals(architectureClass);
    }
}