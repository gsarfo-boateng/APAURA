package org.apaura.util;

/** Central location for APAURA version and citation text. */
public final class Citation {
    public static final String TOOL_NAME = "APAURA";
    public static final String VERSION = "0.10.1";
    public static final String FULL_NAME =
            "Alternative Polyadenylation Architecture, Usage, and Regulatory Analysis";
    public static final String CITATION =
            "Boateng-Sarfo G, Sigor S, Wang L. APAURA: "
                    + "Alternative Polyadenylation Architecture, Usage, and Regulatory Analysis. "
                    + "Version " + VERSION + ", 2026.";

    private Citation() {
        // Utility class.
    }

    public static void printBanner() {
        System.out.println(TOOL_NAME + " v" + VERSION);
        System.out.println(FULL_NAME);
        System.out.println();
    }

    public static void printCitation() {
        System.out.println("Citation:");
        System.out.println(CITATION);
        System.out.println();
    }

    public static void printFooter() {
        System.out.println();
        System.out.println("If you use APAURA, please cite:");
        System.out.println(CITATION);
    }
}
