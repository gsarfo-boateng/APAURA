package org.apaura.plot;

import org.apaura.util.Args;
import org.apaura.util.RunMetadata;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Java entry point for APAURA publication-style locus plots.
 *
 * Java owns command parsing and run metadata. The plot renderer is a bundled R
 * resource so APAURA can produce publication graphics without duplicating the
 * statistical plotting layer in Java.
 */
public final class PlotRegionCommand {

    private static final String RESOURCE =
            "/org/apaura/plot/APAURA_region_plot.R";

    private PlotRegionCommand() {
    }

    public static void run(Args args) throws Exception {
        String resultsArg =
                args.get("results", "");

        boolean standalone =
                isEmpty(resultsArg);

        File results = null;

        Map<String, String> metadata =
                java.util.Collections.emptyMap();

        if (!standalone) {

            results = new File(
                    resultsArg
            ).getCanonicalFile();

            if (!results.isDirectory()) {
                throw new IllegalArgumentException(
                        "Result directory not found: " + results
                );
            }

            metadata =
                    RunMetadata.read(results);
        }


        // ----------------------------------------------------
        // Resolve run information
        // ----------------------------------------------------

        String mode = firstNonEmpty(
                args.get("mode", ""),
                metadata.get("analysis")
        );

        if (isEmpty(mode)) {
            mode = standalone
                    ? "short-read"
                    : detectMode(results);
        }

        mode = normalizeMode(mode);

        if (standalone &&
                !"short-read".equals(mode)) {

            throw new IllegalArgumentException(
                    "Standalone plot-region currently supports "
                            + "short-read BAM/SAM coverage only. "
                            + "Use --results for run-direct/run-isoquant."
            );
        }

        String samples = firstNonEmpty(
                args.get("samples", ""),
                metadata.get("samples")
        );

        String gtf = firstNonEmpty(
                args.get("gtf", ""),
                metadata.get("gtf")
        );

        String control = firstNonEmpty(
                args.get("control", ""),
                metadata.get("control")
        );

        String caseCondition = firstNonEmpty(
                args.get("case", ""),
                metadata.get("case")
        );

        String library = firstNonEmpty(
                args.get("library", ""),
                metadata.get("library"),
                "unstranded"
        );

        String samtools = firstNonEmpty(
                args.get("samtools", ""),
                metadata.get("samtools"),
                "samtools"
        );

        String minMapq = firstNonEmpty(
                args.get("min-mapq", ""),
                metadata.get("min_mapq"),
                "10"
        );

        String rscript = args.get(
                "rscript",
                "Rscript"
        );

        // ----------------------------------------------------
        // Gene(s) / region
        // ----------------------------------------------------

        List<String> genes = parseGenes(
                args.get("gene", ""),
                args.get("genes-file", "")
        );

        String region = args.get("region", "");

        boolean haveGenes = !genes.isEmpty();
        boolean haveRegion = !isEmpty(region);

        if (haveGenes == haveRegion) {
            throw new IllegalArgumentException(
                    "Specify gene(s) with --gene/--genes-file "
                            + "or exactly one --region."
            );
        }

        // ----------------------------------------------------
        // Validate resolved files
        // ----------------------------------------------------

        requireResolved("--samples", samples);
        requireResolved("--gtf", gtf);
        requireResolved("--control", control);
        requireResolved("--case", caseCondition);

        File sampleFile = new File(samples).getCanonicalFile();
        File gtfFile = new File(gtf).getCanonicalFile();

        if (!sampleFile.isFile()) {
            throw new IllegalArgumentException(
                    "Sample sheet not found: " + sampleFile
            );
        }

        if (!gtfFile.isFile()) {
            throw new IllegalArgumentException(
                    "GTF not found: " + gtfFile
            );
        }

        // ----------------------------------------------------
        // Targets
        // ----------------------------------------------------

        List<String> targets = new ArrayList<String>();

        if (haveGenes) {
            targets.addAll(genes);
        } else {
            targets.add(region);
        }

        boolean multi = targets.size() > 1;

        System.err.println(
                "[APAURA] plot-region mode=" + mode
                        + (standalone
                            ? " standalone=true"
                            : "")
                        + " targets=" + targets.size()
                        + " library=" + library
        );

        // ----------------------------------------------------
        // Materialize bundled renderer once
        // ----------------------------------------------------

        File renderer = materializeRenderer();

        try {

            for (int i = 0; i < targets.size(); i++) {

                String target = targets.get(i);

                String gene = haveGenes
                        ? target
                        : "";

                String targetRegion = haveGenes
                        ? ""
                        : region;

                renderTarget(
                        args,
                        results,
                        mode,
                        sampleFile,
                        gtfFile,
                        control,
                        caseCondition,
                        library,
                        samtools,
                        minMapq,
                        rscript,
                        renderer,
                        gene,
                        targetRegion,
                        multi,
                        i + 1,
                        targets.size()
                );
            }

        } finally {

            if (renderer.getName().startsWith(
                    "apaura_plot_region_"
            )) {
                renderer.delete();
            }
        }
    }

    /**
     * Render one gene or one genomic region.
     *
     * Multiple genes are rendered sequentially so the persistent R-side
     * normalization cache can be reused after the first target.
     */
    private static void renderTarget(
            Args args,
            File results,
            String mode,
            File sampleFile,
            File gtfFile,
            String control,
            String caseCondition,
            String library,
            String samtools,
            String minMapq,
            String rscript,
            File renderer,
            String gene,
            String region,
            boolean multi,
            int index,
            int total
    ) throws Exception {

        String target = !isEmpty(gene)
                ? gene
                : region;

        File outFile;

        // ----------------------------------------------------
        // Output
        // ----------------------------------------------------

        if (multi) {

            String outArg = args.get("out", "");

            File outDir = isEmpty(outArg)
                    ? defaultPlotDirectory(results)
                    : new File(outArg);

            outDir = outDir.getCanonicalFile();

            if (!outDir.exists() && !outDir.mkdirs()) {
                throw new IOException(
                        "Cannot create plot output directory: "
                                + outDir
                );
            }

            if (!outDir.isDirectory()) {
                throw new IllegalArgumentException(
                        "For multiple genes, --out must be a directory: "
                                + outDir
                );
            }

            outFile = new File(
                    outDir,
                    sanitize(target) + ".APAURA"
            ).getCanonicalFile();

        } else {

            String defaultOut = new File(
                    defaultPlotDirectory(results),
                    sanitize(target) + ".APAURA"
            ).getPath();

            outFile = new File(
                    args.get("out", defaultOut)
            ).getCanonicalFile();

            File parent = outFile.getParentFile();

            if (parent != null &&
                    !parent.exists() &&
                    !parent.mkdirs()) {

                throw new IOException(
                        "Cannot create plot output directory: "
                                + parent
                );
            }
        }

        // ----------------------------------------------------
        // R command
        // ----------------------------------------------------

        File renderResults = results;

        if (renderResults == null) {
            renderResults =
                    prepareStandaloneResults(outFile);
        }

        List<String> command =
                new ArrayList<String>();

        command.add(rscript);
        command.add(renderer.getAbsolutePath());

        add(command, "--mode", mode);
        add(command, "--results-dir", renderResults.getAbsolutePath());

        if (results == null) {
            command.add("--standalone");
        }
        add(command, "--samples", sampleFile.getCanonicalPath());
        add(command, "--gtf", gtfFile.getCanonicalPath());
        add(command, "--control", control);
        add(command, "--case", caseCondition);
        add(command, "--library", library);
        add(command, "--out", outFile.getPath());

        if (!isEmpty(gene)) {
            add(command, "--gene", gene);
        } else {
            add(command, "--region", region);
        }

        add(
                command,
                "--flank",
                args.get("flank", "1000")
        );

        add(
                command,
                "--bin-size",
                args.get("bin-size", "10")
        );

        add(
                command,
                "--min-mapq",
                minMapq
        );

        add(
                command,
                "--merge",
                args.get("merge", "mean")
        );

        add(
                command,
                "--fdr",
                args.get("fdr", "0.05")
        );

        add(
                command,
                "--apa-effect",
                args.get("min-delta-pau", "0.10")
        );

        add(
                command,
                "--ipa-effect",
                args.get("min-ipa-delta", "0.10")
        );

        add(
                command,
                "--format",
                args.get("format", "all")
        );

        add(
                command,
                "--width",
                args.get("width", "12")
        );

        add(
                command,
                "--height",
                args.get("height", "9")
        );

        add(
                command,
                "--samtools",
                samtools
        );

        add(
                command,
                "--max-transcripts",
                args.get("max-transcripts", "20")
        );

        if (args.flag("show-replicates")) {
            command.add("--show-replicates");
        }
        if (args.flag("show-overview")) {
            command.add("--show-overview");
        }
        if (args.flag("verbose")) {
            command.add("--verbose");
        }

        if (args.flag("no-ipa")) {
            command.add("--no-ipa");
        }

        if (args.flag("keep-duplicates")) {
            command.add("--keep-duplicates");
        }

        if (args.flag("genomic-orientation")) {
            command.add("--genomic-orientation");
        }

        if (args.flag("show-transcripts")) {
            command.add("--show-transcripts");
        }

        // ----------------------------------------------------
        // Run
        // ----------------------------------------------------

        String progress = total > 1
                ? " " + index + "/" + total
                : "";

        System.err.println(
                "[APAURA] plot-region"
                        + progress
                        + " mode=" + mode
                        + " target=" + target
        );

        // BEGIN APAURA_SC_PLOT_REGION_GROUP
        String plotRegionGroup = args.get("group", "");

        if (!isEmpty(plotRegionGroup)) {
            add(
                    command,
                    "--group",
                    plotRegionGroup
            );
        }
        // END APAURA_SC_PLOT_REGION_GROUP

        ProcessBuilder processBuilder =
                new ProcessBuilder(command);

        processBuilder.inheritIO();

        int exit =
                processBuilder.start().waitFor();

        if (exit != 0) {
            throw new IllegalStateException(
                    "APAURA region renderer failed for "
                            + target
                            + " with status "
                            + exit
            );
        }

        System.err.println(
                "[APAURA] plot-region complete: "
                        + outFile
        );
    }

    // ========================================================
    // Gene-list handling
    // ========================================================

    private static File defaultPlotDirectory(
            File results
    ) throws IOException {

        if (results != null) {
            return new File(
                    results,
                    "plots"
            ).getCanonicalFile();
        }

        return new File(
                "plots"
        ).getCanonicalFile();
    }

    /**
     * Standalone plot-region reuses the existing renderer by supplying
     * an empty APA site table plus a persistent normalization cache
     * directory beside the requested plot output.
     */
    private static File prepareStandaloneResults(
            File outFile
    ) throws IOException {

        File parent =
                outFile.getParentFile();

        if (parent == null) {
            parent =
                    new File(".")
                            .getCanonicalFile();
        }

        File cacheDir =
                new File(
                        parent,
                        ".apaura_plot_cache"
                ).getCanonicalFile();

        if (!cacheDir.exists() &&
                !cacheDir.mkdirs()) {

            throw new IOException(
                    "Cannot create standalone plot cache: "
                            + cacheDir
            );
        }

        File apaSites =
                new File(
                        cacheDir,
                        "apa_sites.tsv"
                );

        /*
         * Always rewrite the tiny standalone APA table.
         *
         * A header-only TSV makes readr infer pas_position_1based as
         * character, which breaks dplyr::between() against numeric
         * genomic coordinates. One impossible placeholder row forces
         * numeric inference while being guaranteed to disappear during
         * chromosome filtering.
         */
        BufferedWriter writer =
                new BufferedWriter(
                        new FileWriter(apaSites)
                );

        try {

            writer.write(
                    "pas_position_1based\t"
                            + "delta_pau\t"
                            + "fdr\t"
                            + "status\t"
                            + "gene_id\t"
                            + "gene_name\t"
                            + "chrom"
            );
            writer.newLine();

            writer.write(
                    "0\t"
                            + "0\t"
                            + "1\t"
                            + "standalone_placeholder\t"
                            + ".\t"
                            + ".\t"
                            + "__APAURA_STANDALONE__"
            );
            writer.newLine();

        } finally {
            writer.close();
        }

        return cacheDir;
    }

    private static List<String> parseGenes(
            String geneArg,
            String genesFile
    ) throws IOException {

        LinkedHashSet<String> genes =
                new LinkedHashSet<String>();

        addGenes(
                genes,
                geneArg
        );

        if (!isEmpty(genesFile)) {

            File file =
                    new File(genesFile)
                            .getCanonicalFile();

            if (!file.isFile()) {
                throw new IllegalArgumentException(
                        "Gene list not found: "
                                + file
                );
            }

            BufferedReader reader =
                    new BufferedReader(
                            new FileReader(file)
                    );

            try {

                String line;

                while ((line = reader.readLine()) != null) {

                    // Allow comments in gene files.
                    int comment =
                            line.indexOf('#');

                    if (comment >= 0) {
                        line = line.substring(
                                0,
                                comment
                        );
                    }

                    addGenes(
                            genes,
                            line
                    );
                }

            } finally {
                reader.close();
            }
        }

        return new ArrayList<String>(
                genes
        );
    }

    private static void addGenes(
            LinkedHashSet<String> genes,
            String text
    ) {

        if (isEmpty(text)) {
            return;
        }

        String[] tokens =
                text.split("[,\\s]+");

        for (String token : tokens) {

            String gene =
                    token.trim();

            if (!gene.isEmpty()) {
                genes.add(gene);
            }
        }
    }

    // ========================================================
    // R renderer
    // ========================================================

    private static File materializeRenderer()
            throws IOException {

        InputStream in =
                PlotRegionCommand.class
                        .getResourceAsStream(
                                RESOURCE
                        );

        if (in == null) {

            File fallback =
                    new File(
                            "scripts/APAURA_region_plot.R"
                    );

            if (fallback.isFile()) {
                return fallback.getCanonicalFile();
            }

            throw new IOException(
                    "Bundled APAURA region renderer was not found. "
                            + "Expected classpath resource "
                            + RESOURCE
                            + " or scripts/APAURA_region_plot.R"
            );
        }

        File tmp =
                File.createTempFile(
                        "apaura_plot_region_",
                        ".R"
                );

        tmp.deleteOnExit();

        FileOutputStream out =
                new FileOutputStream(tmp);

        try {

            byte[] buffer =
                    new byte[8192];

            int n;

            while ((n = in.read(buffer)) >= 0) {
                out.write(
                        buffer,
                        0,
                        n
                );
            }

        } finally {

            try {
                in.close();
            } finally {
                out.close();
            }
        }

        return tmp;
    }

    // ========================================================
    // Mode detection
    // ========================================================

    private static String detectMode(
            File results
    ) throws IOException {

        // BEGIN APAURA_SC_PLOT_REGION_MODE
        if (
                new File(
                        results,
                        "apaura_site_usage.tsv"
                ).isFile()
                        &&
                (
                        new File(
                                results,
                                "observed_unreplicated_events.with_architecture.tsv"
                        ).isFile()
                                ||
                        new File(
                                results,
                                "all_results.with_architecture.tsv"
                        ).isFile()
                                ||
                        new File(
                                results,
                                "unified_results.with_architecture.tsv"
                        ).isFile()
                )
        ) {
            return "short-read";
        }
        // END APAURA_SC_PLOT_REGION_MODE

        // Compact short-read output. Prefer the public three-file schema,
        // while retaining compatibility with older short-read result layouts.
        if (
                (
                        new File(results, "apa.tsv").isFile()
                                && new File(results, "apa_events.tsv").isFile()
                                && new File(results, "apa_sites.tsv").isFile()
                )
                        || new File(new File(results, "details"), "pas_coverage.tsv").isFile()
                        || new File(results, "pas_coverage.tsv").isFile()
                        || (new File(results, "pas.tsv").isFile()
                                        && new File(results, "apa_sites.tsv").isFile()
                                        && new File(results, "apa.tsv").isFile())) {
            return "short-read";
        }

        // IsoQuant result metadata.
        File runInfo =
                new File(results, "run_info.txt");

        if (runInfo.isFile()) {

            String text =
                    readSmallFile(runInfo)
                            .toLowerCase(
                                    Locale.ROOT
                            );

            if (text.contains("isoquant")) {
                return "run-isoquant";
            }
        }

        // Direct-RNA output.
        if (
                new File(
                        results,
                        "pas_counts.tsv"
                ).isFile()

                        &&

                        new File(
                                results,
                                "pau.tsv"
                        ).isFile()
        ) {
            return "run-direct";
        }

        throw new IllegalArgumentException(
                "Cannot infer APAURA analysis mode from "
                        + results
                        + ". Supply --mode "
                        + "short-read|run-direct|run-isoquant."
        );
    }

    private static String normalizeMode(
            String mode
    ) {

        if (mode == null) {
            return "";
        }

        String m =
                mode.trim()
                        .toLowerCase(
                                Locale.ROOT
                        );

        if (
                m.equals("short")
                        || m.equals("run-short")
                        || m.equals("short-read")
        ) {
            return "short-read";
        }

        if (
                m.equals("direct")
                        || m.equals("bam")
                        || m.equals("run-direct")
        ) {
            return "run-direct";
        }

        if (
                m.equals("isoquant")
                        || m.equals("run-isoquant")
        ) {
            return "run-isoquant";
        }

        throw new IllegalArgumentException(
                "Unknown plot mode: "
                        + mode
        );
    }

    // ========================================================
    // Helpers
    // ========================================================

    private static String readSmallFile(
            File file
    ) throws IOException {

        BufferedReader reader =
                new BufferedReader(
                        new FileReader(file)
                );

        StringBuilder text =
                new StringBuilder();

        try {

            String line;

            while ((line = reader.readLine()) != null) {

                text.append(line)
                        .append('\n');

                if (text.length() > 200000) {
                    break;
                }
            }

        } finally {
            reader.close();
        }

        return text.toString();
    }

    private static void add(
            List<String> command,
            String key,
            String value
    ) {

        command.add(key);
        command.add(value);
    }

    private static String firstNonEmpty(
            String... values
    ) {

        for (String value : values) {

            if (!isEmpty(value)) {
                return value;
            }
        }

        return "";
    }

    private static boolean isEmpty(
            String value
    ) {

        return value == null
                || value.trim().isEmpty();
    }

    private static void requireResolved(
            String name,
            String value
    ) {

        if (isEmpty(value)) {

            throw new IllegalArgumentException(
                    name
                            + " is required because it was not supplied "
                            + "and is not available in apaura_run.tsv."
            );
        }
    }

    private static String sanitize(
            String value
    ) {

        return value.replaceAll(
                "[^A-Za-z0-9._-]+",
                "_"
        );
    }
}