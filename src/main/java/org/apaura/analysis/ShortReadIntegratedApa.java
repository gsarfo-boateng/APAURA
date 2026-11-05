package org.apaura.analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Integrates short-read TUTR, ALE, IPA and structural architecture at gene level. */
public final class ShortReadIntegratedApa {

    private ShortReadIntegratedApa() {
        // Utility class.
    }

    public static final class Result {
        public ArchitectureCatalog.Entry architecture;

        public int nTutrGroups;
        public int nTutrChangedGroups;
        public String tutrRepresentativeGroup = "NA";
        public boolean tutrChange;
        public String tutrObservedDirection = "not_applicable";
        public String tutrDirection = "not_applicable";
        public double tutrDeltaWpi = Double.NaN;
        public double tutrShiftMagnitude = Double.NaN;
        public double tutrJsd = Double.NaN;
        public double tutrMinPvalue = Double.NaN;
        public double tutrMinFdr = Double.NaN;
        public String tutrStatus = "not_tested";

        public boolean aleChange;
        public String aleObservedDirection = "not_applicable";
        public String aleDirection = "not_applicable";
        public double aleDeltaWpi = Double.NaN;
        public double aleSwitchIndex = Double.NaN;
        public double aleMinPvalue = Double.NaN;
        public double aleMinFdr = Double.NaN;
        public String aleStatus = "not_tested";

        public int nIpaEvents;
        public int nIpaChangedEvents;
        public String ipaRepresentativeEvent = "NA";
        public boolean ipaChange;
        public String ipaDirection = "not_applicable";
        public double ipaDeltaUsage = Double.NaN;
        public double ipaMaxAbsDelta = Double.NaN;
        public double ipaMinPvalue = Double.NaN;
        public double ipaMinFdr = Double.NaN;
        public String ipaStatus = "not_tested";

        public String changedMechanisms = "none";
        public String mechanismDirections = "none";
        public boolean apaChange;
        public String apaDirection = "stable";
        public String apaEvent = "none:stable";
        public double minPvalue = Double.NaN;
        public double minFdr = Double.NaN;
        public String status = "not_tested";
    }

    public static List<Result> integrate(
            List<ArchitectureCatalog.Entry> architecture,
            List<ShortReadDifferential.GroupResult> tutrGroups,
            List<ShortReadAle.GeneResult> aleGenes,
            List<ShortReadIpa.GeneResult> ipaGenes
    ) {
        Map<String, List<ShortReadDifferential.GroupResult>> tutrByGene =
                new LinkedHashMap<String, List<ShortReadDifferential.GroupResult>>();
        for (ShortReadDifferential.GroupResult group : tutrGroups) {
            tutrByGene.computeIfAbsent(
                    group.group.geneId,
                    key -> new ArrayList<ShortReadDifferential.GroupResult>()
            ).add(group);
        }

        Map<String, ShortReadAle.GeneResult> aleByGene =
                new LinkedHashMap<String, ShortReadAle.GeneResult>();
        for (ShortReadAle.GeneResult gene : aleGenes) {
            aleByGene.put(gene.group.geneId, gene);
        }

        Map<String, ShortReadIpa.GeneResult> ipaByGene =
                new LinkedHashMap<String, ShortReadIpa.GeneResult>();
        for (ShortReadIpa.GeneResult gene : ipaGenes) {
            ipaByGene.put(gene.geneId, gene);
        }

        List<Result> output = new ArrayList<Result>();
        for (ArchitectureCatalog.Entry entry : architecture) {
            List<ShortReadDifferential.GroupResult> tutr = tutrByGene.get(entry.geneId);
            ShortReadAle.GeneResult ale = aleByGene.get(entry.geneId);
            ShortReadIpa.GeneResult ipa = ipaByGene.get(entry.geneId);

            boolean terminalArchitecture =
                    ApaArchitecture.isAlternativeTerminalArchitecture(
                            entry.architectureClass
                    );
            if (!terminalArchitecture && ipa == null) {
                continue;
            }

            Result result = new Result();
            result.architecture = entry;
            summarizeTutr(result, tutr);
            summarizeAle(result, ale);
            summarizeIpa(result, ipa);
            finalizeOverall(result);
            output.add(result);
        }

        Collections.sort(output, new Comparator<Result>() {
            @Override
            public int compare(Result left, Result right) {
                int c = left.architecture.chrom.compareTo(right.architecture.chrom);
                return c != 0 ? c : left.architecture.geneId.compareTo(right.architecture.geneId);
            }
        });
        return output;
    }

    private static void summarizeTutr(
            Result result,
            List<ShortReadDifferential.GroupResult> groups
    ) {
        if (groups == null || groups.isEmpty()) {
            return;
        }
        result.nTutrGroups = groups.size();
        ShortReadDifferential.GroupResult representative = null;
        boolean anyOk = false;
        boolean anySparse = false;
        boolean anyEvidenceFailure = false;

        for (ShortReadDifferential.GroupResult group : groups) {
            if (group.apaChange) {
                result.tutrChange = true;
                result.nTutrChangedGroups++;
            }
            result.tutrMinPvalue = finiteMin(result.tutrMinPvalue, group.minPvalue);
            result.tutrMinFdr = finiteMin(result.tutrMinFdr, group.minFdr);
            if (representative == null || betterRepresentative(group, representative)) {
                representative = group;
            }
            if ("ok".equals(group.status)) anyOk = true;
            else if ("insufficient_replicate_support".equals(group.status)) anySparse = true;
            else if ("insufficient_evidence".equals(group.status)) anyEvidenceFailure = true;
        }

        if (representative != null) {
            result.tutrRepresentativeGroup = representative.group.id;
            result.tutrDeltaWpi = representative.deltaWpi;
            result.tutrShiftMagnitude = representative.apd;
            result.tutrJsd = representative.jsd;
            result.tutrObservedDirection = representative.observedDirection;
            result.tutrDirection = representative.apaDirection;
        }
        result.tutrStatus = anyOk
                ? "ok"
                : anySparse
                ? "insufficient_replicate_support"
                : anyEvidenceFailure
                ? "insufficient_evidence"
                : "not_tested";
    }

    private static boolean betterRepresentative(
            ShortReadDifferential.GroupResult candidate,
            ShortReadDifferential.GroupResult current
    ) {
        if (candidate.apaChange != current.apaChange) {
            return candidate.apaChange;
        }

        boolean candidateObserved =
                !"stable".equals(candidate.apaDirection)
                && !"not_applicable".equals(candidate.apaDirection);
        boolean currentObserved =
                !"stable".equals(current.apaDirection)
                && !"not_applicable".equals(current.apaDirection);

        if (candidateObserved != currentObserved) {
            return candidateObserved;
        }

        if (Double.isFinite(candidate.minFdr) && Double.isFinite(current.minFdr)) {
            return candidate.minFdr < current.minFdr;
        }
        if (Double.isFinite(candidate.minFdr)) return true;
        if (Double.isFinite(current.minFdr)) return false;
        double candidateApd = Double.isFinite(candidate.apd) ? candidate.apd : -1.0;
        double currentApd = Double.isFinite(current.apd) ? current.apd : -1.0;
        return candidateApd > currentApd;
    }

    private static void summarizeAle(Result result, ShortReadAle.GeneResult ale) {
        if (ale == null) return;
        result.aleChange = ale.aleChange;
        result.aleObservedDirection = ale.observedDirection;
        result.aleDirection = ale.aleDirection;
        result.aleDeltaWpi = ale.deltaWpi;
        result.aleSwitchIndex = ale.switchIndex;
        result.aleMinPvalue = ale.minPvalue;
        result.aleMinFdr = ale.minFdr;
        result.aleStatus = ale.status;
    }

    private static void summarizeIpa(Result result, ShortReadIpa.GeneResult ipa) {
        if (ipa == null) return;
        result.nIpaEvents = ipa.nIpaEvents;
        result.nIpaChangedEvents = ipa.nSignificantEvents;
        result.ipaRepresentativeEvent = ipa.representativeEvent;
        result.ipaChange = ipa.nSignificantEvents > 0;
        result.ipaDirection = ipa.ipaDirection;
        result.ipaDeltaUsage = ipa.representativeDelta;
        result.ipaMaxAbsDelta = ipa.maxAbsDelta;
        result.ipaMinPvalue = ipa.minPvalue;
        result.ipaMinFdr = ipa.minFdr;
        result.ipaStatus = ipa.status;
    }

    private static void finalizeOverall(Result result) {
        List<String> mechanisms = new ArrayList<String>();
        List<String> directions = new ArrayList<String>();

        // Formal statistical support.
        if (result.ipaChange) {
            mechanisms.add("IPA");
            directions.add("IPA=" + result.ipaDirection);
        }
        if (result.aleChange) {
            mechanisms.add("ALE");
            directions.add("ALE=" + result.aleDirection);
        }
        if (result.tutrChange) {
            mechanisms.add("TUTR");
            directions.add("TUTR=" + result.tutrDirection);
        }

        result.changedMechanisms = join(mechanisms, "+", "none");
        result.mechanismDirections = join(directions, ";", "none");
        result.apaChange = !mechanisms.isEmpty();

        // Observed biological direction is independent of FDR.
        List<String> observedMechanisms = new ArrayList<String>();
        List<String> observedDirections = new ArrayList<String>();

        addObserved(
                observedMechanisms, observedDirections,
                "IPA", result.ipaDirection
        );
        addObserved(
                observedMechanisms, observedDirections,
                "ALE", result.aleDirection
        );
        addObserved(
                observedMechanisms, observedDirections,
                "TUTR", result.tutrDirection
        );

        if (observedDirections.isEmpty()) {
            result.apaDirection = "stable";
        } else if (observedDirections.size() == 1) {
            result.apaDirection = observedDirections.get(0);
        } else {
            result.apaDirection = "complex";
        }

        // apaEvent retains formal-mechanism semantics.
        if (mechanisms.size() == 1) {
            String mechanism = mechanisms.get(0);
            result.apaEvent =
                    mechanism + ":" + mechanismDirection(result, mechanism);
        } else if (mechanisms.size() > 1) {
            result.apaEvent =
                    "mixed".equals(result.architecture.architectureClass)
                            && mechanisms.contains("TUTR")
                            && mechanisms.contains("ALE")
                            ? "mixed:complex"
                            : "multi_mechanism:complex";
        } else {
            result.apaEvent =
                    result.architecture.architectureClass
                            + ":" + result.apaDirection;
        }

        result.minPvalue = finiteMin(
                result.tutrMinPvalue,
                finiteMin(result.aleMinPvalue, result.ipaMinPvalue)
        );
        result.minFdr = finiteMin(
                result.tutrMinFdr,
                finiteMin(result.aleMinFdr, result.ipaMinFdr)
        );

        boolean anyOk = "ok".equals(result.tutrStatus)
                || "ok".equals(result.aleStatus)
                || "ok".equals(result.ipaStatus);
        boolean anySparse =
                "insufficient_replicate_support".equals(result.tutrStatus)
                || "insufficient_replicate_support".equals(result.aleStatus)
                || "insufficient_replicate_support".equals(result.ipaStatus);
        boolean anyEvidenceFailure =
                "insufficient_evidence".equals(result.tutrStatus)
                || "insufficient_evidence".equals(result.aleStatus)
                || "insufficient_evidence".equals(result.ipaStatus);

        result.status = anyOk
                ? "ok"
                : anySparse
                ? "insufficient_replicate_support"
                : anyEvidenceFailure
                ? "insufficient_evidence"
                : "not_tested";
    }

    private static void addObserved(
            List<String> mechanisms,
            List<String> directions,
            String mechanism,
            String direction
    ) {
        if (direction == null
                || "stable".equals(direction)
                || "not_applicable".equals(direction)) {
            return;
        }

        mechanisms.add(mechanism);
        directions.add(direction);
    }

    private static String mechanismDirection(
            Result result,
            String mechanism
    ) {
        if ("TUTR".equals(mechanism)) {
            return result.tutrDirection;
        }
        if ("ALE".equals(mechanism)) {
            return result.aleDirection;
        }
        if ("IPA".equals(mechanism)) {
            return result.ipaDirection;
        }
        return result.apaDirection;
    }

    private static double finiteMin(double left, double right) {
        if (Double.isFinite(left) && Double.isFinite(right)) return Math.min(left, right);
        if (Double.isFinite(left)) return left;
        if (Double.isFinite(right)) return right;
        return Double.NaN;
    }

    private static String join(List<String> values, String separator, String empty) {
        if (values.isEmpty()) return empty;
        StringBuilder out = new StringBuilder();
        for (String value : values) {
            if (out.length() > 0) out.append(separator);
            out.append(value);
        }
        return out.toString();
    }
}
