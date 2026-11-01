package org.apaura.analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/** Registry of supported 3'-end sequencing protocol presets. */
public final class ProtocolRegistry {

    private static final Map<String, ProtocolSpec> SPECS =
            new LinkedHashMap<String, ProtocolSpec>();

    static {
        add("generic",       ProtocolSpec.Boundary.TRANSCRIPT_3P, false, true,  20);
        add("3seq",          ProtocolSpec.Boundary.TRANSCRIPT_3P, false, true,  20);
        add("3p-seq",        ProtocolSpec.Boundary.TRANSCRIPT_3P, false, false, 20);
        add("pas-seq",       ProtocolSpec.Boundary.TRANSCRIPT_3P, false, true,  20);
        add("polyaseq",      ProtocolSpec.Boundary.TRANSCRIPT_3P, false, true,  20);
        add("3reads",        ProtocolSpec.Boundary.TRANSCRIPT_3P, false, false, 20);
        add("3reads+",       ProtocolSpec.Boundary.TRANSCRIPT_3P, false, false, 20);
        add("pat-seq",       ProtocolSpec.Boundary.TRANSCRIPT_3P, false, false, 20);
        add("a-seq",         ProtocolSpec.Boundary.TRANSCRIPT_3P, false, true,  20);
        add("a-seq2",        ProtocolSpec.Boundary.TRANSCRIPT_3P, false, true,  20);
        add("sapas",         ProtocolSpec.Boundary.TRANSCRIPT_3P, false, true,  20);
        add("3t-fill",       ProtocolSpec.Boundary.TRANSCRIPT_3P, false, true,  20);
        add("paperclip",     ProtocolSpec.Boundary.TRANSCRIPT_3P, false, false, 20);
        add("2p-seq",        ProtocolSpec.Boundary.TRANSCRIPT_3P, false, true,  20);
        add("pa-seq",        ProtocolSpec.Boundary.TRANSCRIPT_3P, false, true,  20);
        add("wtts-seq",      ProtocolSpec.Boundary.TRANSCRIPT_3P, false, true,  20);
        add("maps",          ProtocolSpec.Boundary.TRANSCRIPT_3P, false, true,  20);
        add("3pc",           ProtocolSpec.Boundary.TRANSCRIPT_3P, false, true,  20);
        add("exprss",        ProtocolSpec.Boundary.TRANSCRIPT_3P, false, true,  20);
        add("quantseq-fwd",  ProtocolSpec.Boundary.TRANSCRIPT_3P, false, true,  20);
        add("quantseq-rev",  ProtocolSpec.Boundary.ALIGNMENT_5P,  true,  true,  20);
        add("quantseq-ffpe", ProtocolSpec.Boundary.TRANSCRIPT_3P, false, true,  20);
    }

    private ProtocolRegistry() {}

    private static void add(
            String name,
            ProtocolSpec.Boundary boundary,
            boolean flipAlignmentStrand,
            boolean internalPrimingRecommended,
            int defaultClusterDistance
    ) {
        SPECS.put(
                name,
                new ProtocolSpec(
                        name,
                        boundary,
                        flipAlignmentStrand,
                        internalPrimingRecommended,
                        defaultClusterDistance
                )
        );
    }

    public static ProtocolSpec get(String name) {
        String key = normalize(name);
        ProtocolSpec spec = SPECS.get(key);
        if (spec == null) {
            throw new IllegalArgumentException(
                    "Unknown 3'-seq protocol: " + name
                            + ". Supported protocols: "
                            + String.join(", ", names())
            );
        }
        return spec;
    }

    public static List<String> names() {
        return Collections.unmodifiableList(
                new ArrayList<String>(SPECS.keySet())
        );
    }

    private static String normalize(String value) {
        if (value == null) return "";
        String key = value.trim().toLowerCase(Locale.ROOT);
        if ("3pseq".equals(key)) return "3p-seq";
        if ("passeq".equals(key)) return "pas-seq";
        if ("polya-seq".equals(key) || "poly-a-seq".equals(key)) return "polyaseq";
        if ("3readsplus".equals(key) || "3reads-plus".equals(key)) return "3reads+";
        if ("aseq".equals(key)) return "a-seq";
        if ("aseq2".equals(key)) return "a-seq2";
        if ("quantseq".equals(key)) return "quantseq-fwd";
        if ("quantseq-forward".equals(key)) return "quantseq-fwd";
        if ("quantseq-reverse".equals(key)) return "quantseq-rev";
        return key;
    }
}
