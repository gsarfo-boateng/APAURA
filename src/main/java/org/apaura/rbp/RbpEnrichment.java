package org.apaura.rbp;

import org.apaura.stats.FisherExact;
import org.apaura.stats.MultipleTesting;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;

/** Performs architecture-aware RBP motif enrichment around regulated RNA regions. */
public final class RbpEnrichment {

    private RbpEnrichment() {
        // Utility class.
    }

    static final class Gene {
        String id;
        String arch;
        String tutr;
        String ale;
        String ipa;
    }

    static final class Site {
        String id;
        String gene;
        String chrom;
        String region;
        String role;
        int pos0;
        int anchor = -1;
        int rank;
        char strand;
        double delta;
        double p;
        double q;
    }

    static final class Target {
        String id;
        String gene;
        String group;
        String type;
        String chrom;
        String role;
        int pos0;
        int up;
        int down;
        char strand;
        String seq;
        double gc;
        boolean foreground;
    }

    static final class Hit {
        Target target;
        MemeMotif motif;
        int count;
        double best = -1.0;
        final List<Integer> starts = new ArrayList<Integer>();
    }

    static final class Test {
        String group;
        String type;
        String motif;
        String name;
        String backgroundMode;
        int foregroundSites;
        int backgroundSites;
        int foregroundHits;
        int backgroundHits;
        double oddsRatio;
        double pValue;
        double fdr;
        double permutationPValue;
        double permutationFdr;
    }

    public static void run(
            String sitesPath,
            String genesPath,
            String fastaPath,
            String motifPath,
            String outDir,
            int upstream,
            int downstream,
            int spliceUpstream,
            int spliceDownstream,
            double fdr,
            double minDelta,
            double scoreFraction,
            int backgroundRatio,
            int permutations,
            long seed,
            boolean includeNoDirectionalShift
    ) throws Exception {
        Map<String, Gene> genes = readGenes(genesPath);
        List<Site> sites = readSites(sitesPath);
        List<MemeMotif> motifs = MemeParser.read(motifPath);
        File outputDirectory = createOutputDirectory(outDir);

        List<Target> foreground = new ArrayList<Target>();
        List<Target> backgroundCandidates = new ArrayList<Target>();
        Set<String> deduplicationKeys = new HashSet<String>();

        try (FastaReference fasta = new FastaReference(fastaPath)) {
            buildPrimaryTargets(
                    sites,
                    genes,
                    foreground,
                    backgroundCandidates,
                    deduplicationKeys,
                    fasta,
                    upstream,
                    downstream,
                    spliceUpstream,
                    spliceDownstream,
                    fdr,
                    minDelta,
                    includeNoDirectionalShift
            );
        }

        if (foreground.isEmpty()) {
            throw new IllegalArgumentException(
                    "No regulated foreground sites passed --fdr/--min-delta-pau "
                            + "and architecture-direction rules."
            );
        }

        writeTargets(
                new File(outputDirectory, "rbp_windows.tsv"),
                foreground,
                backgroundCandidates
        );

        Map<String, List<Target>> foregroundGroups = groupTargets(foreground);
        List<Hit> allHits = new ArrayList<Hit>();
        List<Test> tests = new ArrayList<Test>();
        Random random = new Random(seed);

        for (Map.Entry<String, List<Target>> entry : foregroundGroups.entrySet()) {
            String[] keyParts = entry.getKey().split("\\|", 2);
            String group = keyParts[0];
            String targetType = keyParts[1];
            List<Target> groupForeground = entry.getValue();

            List<Target> backgroundPool = compatible(
                    backgroundCandidates,
                    group,
                    targetType
            );
            String backgroundMode = "stable_matched";

            if (backgroundPool.isEmpty()) {
                Set<String> foregroundGenes = geneIds(groupForeground);
                backgroundPool = compatibleAllNonForeground(
                        sites,
                        genes,
                        group,
                        targetType,
                        fastaPath,
                        upstream,
                        downstream,
                        spliceUpstream,
                        spliceDownstream,
                        foregroundGenes
                );
                backgroundMode = "nonforeground_fallback";
            }

            List<Target> matchedBackground = matchGc(
                    groupForeground,
                    backgroundPool,
                    backgroundRatio
            );
            if (matchedBackground.isEmpty()) {
                continue;
            }

            testMotifs(
                    group,
                    targetType,
                    groupForeground,
                    matchedBackground,
                    motifs,
                    scoreFraction,
                    permutations,
                    random,
                    backgroundMode,
                    allHits,
                    tests
            );
        }

        adjust(tests);
        writeHits(new File(outputDirectory, "rbp_hits.tsv"), allHits, scoreFraction);
        writeTests(new File(outputDirectory, "rbp.tsv"), tests);
        writePositions(new File(outputDirectory, "rbp_profile.tsv"), allHits);
        writeInfo(
                new File(outputDirectory, "rbp_run_info.txt"),
                motifs.size(),
                foreground.size(),
                backgroundCandidates.size(),
                scoreFraction,
                backgroundRatio,
                permutations,
                seed,
                upstream,
                downstream,
                spliceUpstream,
                spliceDownstream,
                fdr,
                minDelta
        );
    }

    private static File createOutputDirectory(String outDir) throws IOException {
        File directory = new File(outDir);
        if (!directory.exists() && !directory.mkdirs()) {
            throw new IOException("Cannot create " + outDir);
        }
        return directory;
    }

    private static void buildPrimaryTargets(
            List<Site> sites,
            Map<String, Gene> genes,
            List<Target> foreground,
            List<Target> backgroundCandidates,
            Set<String> deduplicationKeys,
            FastaReference fasta,
            int upstream,
            int downstream,
            int spliceUpstream,
            int spliceDownstream,
            double fdr,
            double minDelta,
            boolean includeNoDirectionalShift
    ) throws Exception {
        for (Site site : sites) {
            Gene gene = genes.get(site.gene);
            if (gene == null) {
                continue;
            }

            boolean significant = site.q <= fdr && Math.abs(site.delta) >= minDelta;
            addTutrTargets(
                    site,
                    gene,
                    significant,
                    includeNoDirectionalShift,
                    foreground,
                    backgroundCandidates,
                    deduplicationKeys,
                    fasta,
                    upstream,
                    downstream
            );
            addAleTargets(
                    site,
                    gene,
                    significant,
                    foreground,
                    backgroundCandidates,
                    deduplicationKeys,
                    fasta,
                    upstream,
                    downstream,
                    spliceUpstream,
                    spliceDownstream
            );
            addIpaTargets(
                    site,
                    gene,
                    significant,
                    foreground,
                    backgroundCandidates,
                    deduplicationKeys,
                    fasta,
                    upstream,
                    downstream
            );
        }
    }

    private static void addTutrTargets(
            Site site,
            Gene gene,
            boolean significant,
            boolean includeNoDirectionalShift,
            List<Target> foreground,
            List<Target> backgroundCandidates,
            Set<String> deduplicationKeys,
            FastaReference fasta,
            int upstream,
            int downstream
    ) throws Exception {
        if (!site.role.contains("TUTR")) {
            return;
        }

        String direction = gene.tutr;
        if (validTutr(direction, includeNoDirectionalShift)) {
            if (significant && site.delta > 0) {
                addTarget(
                        foreground,
                        deduplicationKeys,
                        site,
                        "TUTR:" + direction,
                        "PAS",
                        site.pos0,
                        upstream,
                        downstream,
                        true,
                        fasta
                );
            }
        } else if (isStable(direction)) {
            addTarget(
                    backgroundCandidates,
                    deduplicationKeys,
                    site,
                    "TUTR:bg",
                    "PAS",
                    site.pos0,
                    upstream,
                    downstream,
                    false,
                    fasta
            );
        }
    }

    private static void addAleTargets(
            Site site,
            Gene gene,
            boolean significant,
            List<Target> foreground,
            List<Target> backgroundCandidates,
            Set<String> deduplicationKeys,
            FastaReference fasta,
            int upstream,
            int downstream,
            int spliceUpstream,
            int spliceDownstream
    ) throws Exception {
        if (!site.role.contains("ALE")) {
            return;
        }

        String direction = gene.ale;
        if (validAle(direction)) {
            if (!significant || site.delta <= 0) {
                return;
            }
            addTarget(
                    foreground,
                    deduplicationKeys,
                    site,
                    "ALE:" + direction,
                    "PAS",
                    site.pos0,
                    upstream,
                    downstream,
                    true,
                    fasta
            );
            if (site.anchor >= 0) {
                addTarget(
                        foreground,
                        deduplicationKeys,
                        site,
                        "ALE:" + direction,
                        "ALE_acceptor",
                        site.anchor,
                        spliceUpstream,
                        spliceDownstream,
                        true,
                        fasta
                );
            }
        } else if (isStable(direction)) {
            addTarget(
                    backgroundCandidates,
                    deduplicationKeys,
                    site,
                    "ALE:bg",
                    "PAS",
                    site.pos0,
                    upstream,
                    downstream,
                    false,
                    fasta
            );
            if (site.anchor >= 0) {
                addTarget(
                        backgroundCandidates,
                        deduplicationKeys,
                        site,
                        "ALE:bg",
                        "ALE_acceptor",
                        site.anchor,
                        spliceUpstream,
                        spliceDownstream,
                        false,
                        fasta
                );
            }
        }
    }

    private static void addIpaTargets(
            Site site,
            Gene gene,
            boolean significant,
            List<Target> foreground,
            List<Target> backgroundCandidates,
            Set<String> deduplicationKeys,
            FastaReference fasta,
            int upstream,
            int downstream
    ) throws Exception {
        if (!"IPA".equals(site.role)) {
            return;
        }

        String direction = gene.ipa;
        if ("gain".equals(direction) || "loss".equals(direction)) {
            boolean desiredSign = "gain".equals(direction)
                    ? site.delta > 0
                    : site.delta < 0;
            if (significant && desiredSign) {
                addTarget(
                        foreground,
                        deduplicationKeys,
                        site,
                        "IPA:" + direction,
                        "PAS",
                        site.pos0,
                        upstream,
                        downstream,
                        true,
                        fasta
                );
            }
        } else if (isStable(direction)) {
            addTarget(
                    backgroundCandidates,
                    deduplicationKeys,
                    site,
                    "IPA:bg",
                    "PAS",
                    site.pos0,
                    upstream,
                    downstream,
                    false,
                    fasta
            );
        }
    }

    private static void testMotifs(
            String group,
            String targetType,
            List<Target> foreground,
            List<Target> background,
            List<MemeMotif> motifs,
            double scoreFraction,
            int permutations,
            Random random,
            String backgroundMode,
            List<Hit> allHits,
            List<Test> tests
    ) {
        for (MemeMotif motif : motifs) {
            int foregroundHits = scanTargets(
                    foreground,
                    motif,
                    scoreFraction,
                    allHits
            );
            int backgroundHits = scanTargets(
                    background,
                    motif,
                    scoreFraction,
                    allHits
            );

            Test test = new Test();
            test.group = group;
            test.type = targetType;
            test.motif = motif.id;
            test.name = motif.name;
            test.backgroundMode = backgroundMode;
            test.foregroundSites = foreground.size();
            test.backgroundSites = background.size();
            test.foregroundHits = foregroundHits;
            test.backgroundHits = backgroundHits;
            test.oddsRatio = oddsRatio(
                    foregroundHits,
                    foreground.size(),
                    backgroundHits,
                    background.size()
            );
            test.pValue = FisherExact.rightTail(
                    foregroundHits,
                    foreground.size() - foregroundHits,
                    backgroundHits,
                    background.size() - backgroundHits
            );
            test.permutationPValue = permutationPValue(
                    foreground,
                    background,
                    motif,
                    scoreFraction,
                    permutations,
                    random
            );
            tests.add(test);
        }
    }

    private static int scanTargets(
            List<Target> targets,
            MemeMotif motif,
            double scoreFraction,
            List<Hit> allHits
    ) {
        int hitSites = 0;
        for (Target target : targets) {
            Hit hit = scan(target, motif, scoreFraction);
            allHits.add(hit);
            if (hit.count > 0) {
                hitSites++;
            }
        }
        return hitSites;
    }

    private static double oddsRatio(
            int foregroundHits,
            int foregroundSize,
            int backgroundHits,
            int backgroundSize
    ) {
        return ((foregroundHits + 0.5)
                * (backgroundSize - backgroundHits + 0.5))
                / ((foregroundSize - foregroundHits + 0.5)
                * (backgroundHits + 0.5));
    }

    private static Set<String> geneIds(List<Target> targets) {
        Set<String> genes = new HashSet<String>();
        for (Target target : targets) {
            genes.add(target.gene);
        }
        return genes;
    }

    private static List<Target> compatibleAllNonForeground(
            List<Site> sites,
            Map<String, Gene> genes,
            String group,
            String type,
            String fastaPath,
            int upstream,
            int downstream,
            int spliceUpstream,
            int spliceDownstream,
            Set<String> skipGenes
    ) throws Exception {
        List<Target> targets = new ArrayList<Target>();
        Set<String> deduplicationKeys = new HashSet<String>();
        String mechanism = group.split(":", 2)[0];

        try (FastaReference fasta = new FastaReference(fastaPath)) {
            for (Site site : sites) {
                Gene gene = genes.get(site.gene);
                if (gene == null || skipGenes.contains(site.gene)) {
                    continue;
                }

                if ("TUTR".equals(mechanism) && site.role.contains("TUTR")) {
                    addTarget(
                            targets,
                            deduplicationKeys,
                            site,
                            "TUTR:bg",
                            type,
                            site.pos0,
                            upstream,
                            downstream,
                            false,
                            fasta
                    );
                } else if ("ALE".equals(mechanism) && site.role.contains("ALE")) {
                    if ("PAS".equals(type)) {
                        addTarget(
                                targets,
                                deduplicationKeys,
                                site,
                                "ALE:bg",
                                type,
                                site.pos0,
                                upstream,
                                downstream,
                                false,
                                fasta
                        );
                    } else if (site.anchor >= 0) {
                        addTarget(
                                targets,
                                deduplicationKeys,
                                site,
                                "ALE:bg",
                                type,
                                site.anchor,
                                spliceUpstream,
                                spliceDownstream,
                                false,
                                fasta
                        );
                    }
                } else if ("IPA".equals(mechanism) && "IPA".equals(site.role)) {
                    addTarget(
                            targets,
                            deduplicationKeys,
                            site,
                            "IPA:bg",
                            type,
                            site.pos0,
                            upstream,
                            downstream,
                            false,
                            fasta
                    );
                }
            }
        }
        return targets;
    }

    private static boolean validTutr(String direction, boolean includeNoDirectionalShift) {
        return "shortening".equals(direction)
                || "lengthening".equals(direction)
                || (includeNoDirectionalShift
                && "redistribution".equals(direction));
    }

    private static boolean validAle(String direction) {
        return "shortening".equals(direction)
                || "lengthening".equals(direction)
                || "switch".equals(direction);
    }

    private static boolean isStable(String direction) {
        return direction == null
                || "stable".equals(direction)
                || "not_applicable".equals(direction);
    }

    private static void addTarget(
            List<Target> output,
            Set<String> deduplicationKeys,
            Site site,
            String group,
            String type,
            int position0,
            int upstream,
            int downstream,
            boolean foreground,
            FastaReference fasta
    ) throws Exception {
        String key = (foreground ? "F" : "B")
                + "|" + group
                + "|" + type
                + "|" + site.chrom
                + "|" + position0
                + "|" + site.strand
                + "|" + site.gene;
        if (!deduplicationKeys.add(key)) {
            return;
        }

        Target target = new Target();
        target.id = site.id;
        target.gene = site.gene;
        target.group = group;
        target.type = type;
        target.chrom = site.chrom;
        target.pos0 = position0;
        target.strand = site.strand;
        target.role = site.role;
        target.up = upstream;
        target.down = downstream;
        target.foreground = foreground;
        target.seq = fasta.window(
                target.chrom,
                target.pos0,
                target.strand,
                upstream,
                downstream
        );
        target.gc = gcFraction(target.seq);
        output.add(target);
    }

    private static double gcFraction(String sequence) {
        int validBases = 0;
        int gcBases = 0;
        for (char base : sequence.toCharArray()) {
            if (base == 'A' || base == 'C' || base == 'G' || base == 'T') {
                validBases++;
                if (base == 'C' || base == 'G') {
                    gcBases++;
                }
            }
        }
        return validBases == 0 ? 0.0 : (double) gcBases / validBases;
    }

    private static Map<String, List<Target>> groupTargets(List<Target> targets) {
        Map<String, List<Target>> groups = new LinkedHashMap<String, List<Target>>();
        for (Target target : targets) {
            String key = target.group + "|" + target.type;
            groups.computeIfAbsent(
                    key,
                    ignored -> new ArrayList<Target>()
            ).add(target);
        }
        return groups;
    }

    private static List<Target> compatible(
            List<Target> targets,
            String group,
            String type
    ) {
        String mechanism = group.split(":", 2)[0];
        List<Target> compatible = new ArrayList<Target>();
        for (Target target : targets) {
            if (target.type.equals(type)
                    && target.group.startsWith(mechanism + ":")) {
                compatible.add(target);
            }
        }
        return compatible;
    }

    private static List<Target> matchGc(
            List<Target> foreground,
            List<Target> pool,
            int ratio
    ) {
        List<Target> matched = new ArrayList<Target>();
        Set<String> used = new HashSet<String>();

        for (Target foregroundTarget : foreground) {
            List<Target> candidates = new ArrayList<Target>();
            for (Target backgroundTarget : pool) {
                if (!backgroundTarget.gene.equals(foregroundTarget.gene)) {
                    candidates.add(backgroundTarget);
                }
            }

            candidates.sort(Comparator.comparingDouble(
                    candidate -> Math.abs(candidate.gc - foregroundTarget.gc)
            ));

            int selected = 0;
            for (Target candidate : candidates) {
                String key = candidate.type
                        + "|" + candidate.chrom
                        + "|" + candidate.pos0
                        + "|" + candidate.strand
                        + "|" + candidate.gene;
                if (used.add(key)) {
                    matched.add(candidate);
                    selected++;
                    if (selected >= ratio) {
                        break;
                    }
                }
            }
        }
        return matched;
    }

    private static Hit scan(Target target, MemeMotif motif, double scoreFraction) {
        Hit hit = new Hit();
        hit.target = target;
        hit.motif = motif;

        for (int offset = 0;
                offset + motif.width <= target.seq.length();
                offset++) {
            double score = motif.normalizedScore(target.seq, offset);
            if (score >= scoreFraction) {
                hit.count++;
                hit.best = Math.max(hit.best, score);
                hit.starts.add(offset - target.up);
            }
        }
        return hit;
    }

    private static double permutationPValue(
            List<Target> foreground,
            List<Target> background,
            MemeMotif motif,
            double scoreFraction,
            int permutations,
            Random random
    ) {
        if (permutations <= 0) {
            return Double.NaN;
        }

        boolean[] hits = motifPresence(
                foreground,
                background,
                motif,
                scoreFraction
        );
        int foregroundHits = countTrue(hits, 0, foreground.size());
        int backgroundHits = countTrue(hits, foreground.size(), hits.length);
        double observed = enrichmentDifference(
                foregroundHits,
                foreground.size(),
                backgroundHits,
                background.size()
        );

        int[] indices = new int[hits.length];
        for (int index = 0; index < indices.length; index++) {
            indices[index] = index;
        }

        int greaterOrEqual = 0;
        for (int permutation = 0; permutation < permutations; permutation++) {
            shuffle(indices, random);
            int permutedForegroundHits = 0;
            int permutedBackgroundHits = 0;

            for (int index = 0; index < foreground.size(); index++) {
                if (hits[indices[index]]) {
                    permutedForegroundHits++;
                }
            }
            for (int index = foreground.size(); index < indices.length; index++) {
                if (hits[indices[index]]) {
                    permutedBackgroundHits++;
                }
            }

            double difference = enrichmentDifference(
                    permutedForegroundHits,
                    foreground.size(),
                    permutedBackgroundHits,
                    background.size()
            );
            if (difference >= observed - 1.0E-12) {
                greaterOrEqual++;
            }
        }
        return (greaterOrEqual + 1.0) / (permutations + 1.0);
    }

    private static boolean[] motifPresence(
            List<Target> foreground,
            List<Target> background,
            MemeMotif motif,
            double scoreFraction
    ) {
        boolean[] hits = new boolean[foreground.size() + background.size()];
        int index = 0;
        for (Target target : foreground) {
            hits[index++] = scan(target, motif, scoreFraction).count > 0;
        }
        for (Target target : background) {
            hits[index++] = scan(target, motif, scoreFraction).count > 0;
        }
        return hits;
    }

    private static int countTrue(boolean[] values, int start, int end) {
        int count = 0;
        for (int index = start; index < end; index++) {
            if (values[index]) {
                count++;
            }
        }
        return count;
    }

    private static double enrichmentDifference(
            int foregroundHits,
            int foregroundSize,
            int backgroundHits,
            int backgroundSize
    ) {
        return (double) foregroundHits / foregroundSize
                - (double) backgroundHits / backgroundSize;
    }

    private static void shuffle(int[] values, Random random) {
        for (int index = values.length - 1; index > 0; index--) {
            int swapIndex = random.nextInt(index + 1);
            int temporary = values[index];
            values[index] = values[swapIndex];
            values[swapIndex] = temporary;
        }
    }

    private static void adjust(List<Test> tests) {
        double[] fisherPValues = new double[tests.size()];
        double[] permutationPValues = new double[tests.size()];

        for (int index = 0; index < tests.size(); index++) {
            Test test = tests.get(index);
            fisherPValues[index] = test.pValue;
            permutationPValues[index] = Double.isNaN(test.permutationPValue)
                    ? 1.0
                    : test.permutationPValue;
        }

        double[] fisherFdr = MultipleTesting.bh(fisherPValues);
        double[] permutationFdr = MultipleTesting.bh(permutationPValues);

        for (int index = 0; index < tests.size(); index++) {
            Test test = tests.get(index);
            test.fdr = fisherFdr[index];
            test.permutationFdr = Double.isNaN(test.permutationPValue)
                    ? Double.NaN
                    : permutationFdr[index];
        }
    }

    private static File resolveRbpInput(
            String path,
            String... requiredColumns
    ) throws IOException {
        File requested = new File(path);

        if (hasColumns(requested, requiredColumns)) {
            return requested;
        }

        File parent = requested.getAbsoluteFile().getParentFile();
        File detailed = parent == null
                ? null
                : new File(new File(parent, "details"), requested.getName());

        if (detailed != null
                && detailed.isFile()
                && hasColumns(detailed, requiredColumns)) {
            return detailed;
        }

        throw new IllegalArgumentException(
                "RBP enrichment requires detailed APAURA metadata. "
                        + "Expected the required columns in "
                        + requested.getPath()
                        + " or "
                        + (detailed == null ? "details/" : detailed.getPath())
        );
    }

    private static boolean hasColumns(
            File file,
            String... requiredColumns
    ) throws IOException {
        if (!file.isFile()) {
            return false;
        }

        try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
            String header = reader.readLine();
            if (header == null) {
                return false;
            }

            Map<String, Integer> columns = index(header);
            for (String column : requiredColumns) {
                if (!columns.containsKey(column)) {
                    return false;
                }
            }
        }

        return true;
    }

    private static Map<String, Gene> readGenes(String path) throws IOException {
        Map<String, Gene> genes = new HashMap<String, Gene>();
        File source = resolveRbpInput(
                path,
                "gene_id",
                "tutr_direction",
                "ale_direction",
                "ipa_direction"
        );

        try (BufferedReader reader = new BufferedReader(new FileReader(source))) {
            Map<String, Integer> columns = index(reader.readLine());
            String line;
            while ((line = reader.readLine()) != null) {
                String[] fields = line.split("\t", -1);
                Gene gene = new Gene();
                gene.id = value(fields, columns, "gene_id");
                gene.arch = value(fields, columns, "apa_architecture");
                gene.tutr = value(fields, columns, "tutr_direction");
                gene.ale = value(fields, columns, "ale_direction");
                gene.ipa = value(fields, columns, "ipa_direction");
                genes.put(gene.id, gene);
            }
        }
        return genes;
    }

    private static List<Site> readSites(String path) throws IOException {
        List<Site> sites = new ArrayList<Site>();
        File source = resolveRbpInput(
                path,
                "site_id",
                "pos0",
                "architecture_role",
                "terminal_exon_anchor0",
                "site_rank",
                "delta_pau"
        );

        try (BufferedReader reader = new BufferedReader(new FileReader(source))) {
            Map<String, Integer> columns = index(reader.readLine());
            String line;
            while ((line = reader.readLine()) != null) {
                String[] fields = line.split("\t", -1);
                Site site = new Site();
                site.id = value(fields, columns, "site_id");
                site.gene = value(fields, columns, "gene_id");
                site.chrom = value(fields, columns, "chrom");
                site.pos0 = Integer.parseInt(value(fields, columns, "pos0"));
                site.strand = value(fields, columns, "strand").charAt(0);
                site.region = value(fields, columns, "region");
                site.role = value(fields, columns, "architecture_role");

                String anchor = value(fields, columns, "terminal_exon_anchor0");
                if (!"NA".equals(anchor) && !anchor.isEmpty()) {
                    site.anchor = Integer.parseInt(anchor);
                }

                site.rank = Integer.parseInt(value(fields, columns, "site_rank"));
                site.delta = Double.parseDouble(value(fields, columns, "delta_pau"));
                site.p = Double.parseDouble(value(fields, columns, "pvalue"));
                site.q = Double.parseDouble(value(fields, columns, "fdr"));
                sites.add(site);
            }
        }
        return sites;
    }

    private static Map<String, Integer> index(String header) {
        String[] fields = header.split("\t");
        Map<String, Integer> columns = new HashMap<String, Integer>();
        for (int index = 0; index < fields.length; index++) {
            columns.put(fields[index], index);
        }
        return columns;
    }

    private static String value(
            String[] fields,
            Map<String, Integer> columns,
            String name
    ) {
        Integer index = columns.get(name);
        if (index == null) {
            throw new IllegalArgumentException("Missing column " + name);
        }
        return fields[index];
    }

    private static void writeTargets(
            File file,
            List<Target> foreground,
            List<Target> background
    ) throws IOException {
        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println(
                    "set\tgroup\ttarget_type\ttarget_id\tgene_id\tchrom\t"
                            + "anchor_pos0\tstrand\tgc\tsequence"
            );
            for (Target target : foreground) {
                writeTarget(writer, target);
            }
            for (Target target : background) {
                writeTarget(writer, target);
            }
        }
    }

    private static void writeTarget(PrintWriter writer, Target target) {
        writer.printf(
                Locale.ROOT,
                "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%c\t%.4f\t%s%n",
                target.foreground ? "foreground" : "background_candidate",
                target.group,
                target.type,
                target.id,
                target.gene,
                target.chrom,
                target.pos0,
                target.strand,
                target.gc,
                target.seq
        );
    }

    private static void writeHits(
            File file,
            List<Hit> hits,
            double threshold
    ) throws IOException {
        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println(
                    "group\ttarget_type\ttarget_id\tgene_id\tmotif_id\trbp_name\t"
                            + "n_hits\tbest_normalized_score\trelative_hit_starts\tthreshold"
            );
            for (Hit hit : hits) {
                writer.printf(
                        Locale.ROOT,
                        "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%.3f%n",
                        hit.target.group,
                        hit.target.type,
                        hit.target.id,
                        hit.target.gene,
                        hit.motif.id,
                        hit.motif.name,
                        hit.count,
                        hit.best < 0
                                ? "NA"
                                : String.format(Locale.ROOT, "%.4f", hit.best),
                        join(hit.starts),
                        threshold
                );
            }
        }
    }

    private static String join(List<Integer> values) {
        StringBuilder joined = new StringBuilder();
        for (int index = 0; index < values.size(); index++) {
            if (index > 0) {
                joined.append(',');
            }
            joined.append(values.get(index));
        }
        return joined.length() == 0 ? "NA" : joined.toString();
    }

    private static void writeTests(File file, List<Test> tests) throws IOException {
        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println(
                    "group\ttarget_type\tmotif_id\trbp_name\tforeground_sites\t"
                            + "background_sites\tforeground_hit_sites\tbackground_hit_sites\t"
                            + "foreground_hit_rate\tbackground_hit_rate\todds_ratio\t"
                            + "fisher_pvalue\tfisher_fdr\tpermutation_pvalue\t"
                            + "permutation_fdr\tbackground_mode"
            );
            for (Test test : tests) {
                writer.printf(
                        Locale.ROOT,
                        "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.6f\t%.6f\t%.6g\t"
                                + "%.8g\t%.8g\t%s\t%s\t%s%n",
                        test.group,
                        test.type,
                        test.motif,
                        test.name,
                        test.foregroundSites,
                        test.backgroundSites,
                        test.foregroundHits,
                        test.backgroundHits,
                        (double) test.foregroundHits / test.foregroundSites,
                        (double) test.backgroundHits / test.backgroundSites,
                        test.oddsRatio,
                        test.pValue,
                        test.fdr,
                        Double.isNaN(test.permutationPValue)
                                ? "NA"
                                : String.format(
                                        Locale.ROOT,
                                        "%.8g",
                                        test.permutationPValue
                                ),
                        Double.isNaN(test.permutationFdr)
                                ? "NA"
                                : String.format(
                                        Locale.ROOT,
                                        "%.8g",
                                        test.permutationFdr
                                ),
                        test.backgroundMode
                );
            }
        }
    }

    private static void writePositions(File file, List<Hit> hits) throws IOException {
        Map<String, Integer> counts = new TreeMap<String, Integer>();
        for (Hit hit : hits) {
            for (Integer position : hit.starts) {
                String key = hit.target.group
                        + "\t" + hit.target.type
                        + "\t" + hit.motif.id
                        + "\t" + hit.motif.name
                        + "\t" + position;
                Integer previous = counts.get(key);
                counts.put(key, previous == null ? 1 : previous + 1);
            }
        }

        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println(
                    "group\ttarget_type\tmotif_id\trbp_name\trelative_start\thit_count"
            );
            for (Map.Entry<String, Integer> entry : counts.entrySet()) {
                writer.println(entry.getKey() + "\t" + entry.getValue());
            }
        }
    }

    private static void writeInfo(
            File file,
            int motifCount,
            int foregroundCount,
            int backgroundCandidateCount,
            double threshold,
            int ratio,
            int permutations,
            long seed,
            int upstream,
            int downstream,
            int spliceUpstream,
            int spliceDownstream,
            double fdr,
            double delta
    ) throws IOException {
        try (PrintWriter writer = new PrintWriter(file)) {
            writer.println("APAURA RBP motif enrichment");
            writer.println("motifs=" + motifCount);
            writer.println("foreground_targets=" + foregroundCount);
            writer.println("stable_background_candidates=" + backgroundCandidateCount);
            writer.println("normalized_motif_score_threshold=" + threshold);
            writer.println("background_ratio=" + ratio);
            writer.println("permutations=" + permutations);
            writer.println("seed=" + seed);
            writer.println("pas_window=-" + upstream + "..+" + downstream);
            writer.println(
                    "ale_acceptor_window=-" + spliceUpstream + "..+" + spliceDownstream
            );
            writer.println("site_fdr=" + fdr);
            writer.println("min_abs_delta_pau=" + delta);
            writer.println(
                    "motif_scanner=independent PWM log2-odds scanner; RNA U is treated as T"
            );
            writer.println(
                    "enrichment_test=one-sided Fisher exact test plus "
                            + "label-permutation enrichment test"
            );
        }
    }
}
