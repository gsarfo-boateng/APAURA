package org.apaura;

import org.apaura.analysis.DirectApaPipeline;
import org.apaura.analysis.EndExtractor;
import org.apaura.analysis.ShortReadPipeline;
import org.apaura.model.EndKey;
import org.apaura.rbp.RbpEnrichment;
import org.apaura.util.Args;
import org.apaura.util.Citation;
import org.apaura.plot.PlotRegionCommand;
import org.apaura.util.RunMetadata;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import org.apaura.sc.Sc3utrCommand;
import java.util.Arrays;

public final class APAURA {
    public static final String VERSION = Citation.VERSION;

    private APAURA() {
    }

    public static void main(String[] argv) {
        try {
            Citation.printBanner();

            if (argv.length == 0 || isHelp(argv[0])) {
                printHelp();
                Citation.printFooter();
                return;
            }

            String command = argv[0];

            if (argv.length > 1 && isHelp(argv[1])) {
                printCommandHelp(command);
                Citation.printFooter();
                return;
            }
            /*
             * Single-cell commands receive their original raw arguments.
             * Do not pass them through APAURA's Args parser.
             */

            if ("sc".equals(command)) {
                runSc(
                        Arrays.copyOfRange(
                                argv,
                                1,
                                argv.length
                        )
                );
                Citation.printFooter();
                return;
            }

            Args args = new Args(argv, 1);
            runCommand(command, args);
            Citation.printFooter();

        } catch (Exception exception) {
            System.err.println("ERROR: " + exception.getMessage());

            if (System.getenv("APAURA_DEBUG") != null) {
                exception.printStackTrace();
            }

            System.exit(2);
        }
    }

    private static void runSc(String[] args) throws Exception {

        if (args.length == 0 || isHelp(args[0])) {
            printScHelp();
            return;
        }

        String subcommand = args[0];

        if ("gtf-sites".equals(subcommand)) {

            if (args.length > 1 && isHelp(args[1])) {
                printScGtfSitesHelp();
                return;
            }

            Sc3utrCommand.runGtfSites(args);
            return;
        }

        if ("import-sites".equals(subcommand)) {

            if (args.length > 1 && isHelp(args[1])) {
                printScImportSitesHelp();
                return;
            }

            Sc3utrCommand.runCore(args);
            return;
        }

        if ("merge-sites".equals(subcommand)) {

            if (args.length > 1 && isHelp(args[1])) {
                printScMergeSitesHelp();
                return;
            }

            Sc3utrCommand.runCore(args);
            return;
        }
        if ("quantify".equals(subcommand)) {

            if (args.length > 1 && isHelp(args[1])) {
                printScQuantifyHelp();
                return;
            }

            Sc3utrCommand.runCore(args);
            return;
        }

        if ("reference".equals(subcommand)) {

            if (args.length > 1 && isHelp(args[1])) {
                printScReferenceHelp();
                return;
            }

            Sc3utrCommand.runReference(args);
            return;
        }

        if ("compare".equals(subcommand)) {

            if (args.length > 1 && isHelp(args[1])) {
                printScCompareHelp();
                return;
            }

            Sc3utrCommand.runCompare(args);
            return;
        }
        throw new IllegalArgumentException(
                "Unknown sc command: " + subcommand
        );
    }

    private static void runCommand(String command, Args args) throws Exception {
        switch (command) {
            case "version":
            case "--version":
            case "-v":
                System.out.println("APAURA " + VERSION);
                return;

            case "citation":
            case "--citation":
                Citation.printCitation();
                return;

            case "ends":
                runEnds(args);
                return;
	    case "plot-region":
                PlotRegionCommand.run(args);
                return;
            case "run-direct":
                runDirect(args, "bam");
                return;

            case "run-isoquant":
                runDirect(args, "isoquant");
                return;

            case "short-read":
            case "run-short":
                runShortRead(args);
                return;

            case "rbp-enrich":
                runRbpEnrichment(args);
                return;

            default:
                throw new IllegalArgumentException("Unknown command: " + command);
        }
    }

    private static boolean isHelp(String value) {
        return value.equals("help") || value.equals("--help") || value.equals("-h");
    }

    private static void runEnds(Args args) throws Exception {
        String alignmentPath = args.require("bam");
        String outputPath = args.require("out");

        Map<EndKey, Integer> counts = EndExtractor.extract(
                alignmentPath,
                args.getInt("min-mapq", 10),
                args.flag("keep-duplicates"),
                args.get("strand-mode", "alignment")
        );

        List<EndKey> ends = new ArrayList<EndKey>(counts.keySet());
        Collections.sort(ends);

        try (PrintWriter writer = new PrintWriter(
                new BufferedWriter(new FileWriter(outputPath)))) {
            writer.println("#chrom\tstart\tend\tstrand\tcount");
            for (EndKey end : ends) {
                writer.println(
                        end.chrom
                                + "\t" + end.pos0
                                + "\t" + (end.pos0 + 1)
                                + "\t" + end.strand.symbol()
                                + "\t" + counts.get(end)
                );
            }
        }
    }

    private static void runDirect(Args args, String defaultInputFormat) throws Exception {
        DirectApaPipeline.run(
                args.require("samples"),
                args.require("gtf"),
                args.require("out"),
                args.get("control", "control"),
                args.get("case", "treated"),
                args.getInt("min-mapq", 10),
                args.flag("keep-duplicates"),
                args.get("strand-mode", "alignment"),
                args.getInt("cluster-distance", 20),
                args.getInt("annotation-slop", 1000),
                args.getInt("min-total", 1),
                args.getDouble("min-delta-pau", 0.10),
                args.getDouble("fdr", 0.05),
                args.getDouble("wpi-threshold", 0.10),
                args.get("input-format", defaultInputFormat),
                args.get("cluster-representative", "mode"),
                args.getInt("min-gene-reads", 10),
                args.getInt("min-site-reads", 5),
                args.getDouble("min-pau", 0.01),
                args.getDouble(
                        "min-internal-exon-delta",
                        args.getDouble("min-architecture-shift", 0.10)
                ),
                args.getDouble("min-ale-switch", 0.10),
                args.getDouble("min-ipa-delta", 0.10),
                args.getDouble("min-tutr-shift", 0.10)
        );
    }

    private static void runShortRead(Args args) throws Exception {
        String mode = args.get("mode", "annotated");

        if (!mode.equalsIgnoreCase("annotated")) {
            throw new IllegalArgumentException(
                    "Only --mode annotated is enabled for conventional RNA-seq in this release."
            );
        }

        ShortReadPipeline.run(
                args.require("samples"),
                args.require("gtf"),
                args.require("out"),
                args.get("control", "control"),
                args.get("case", "treated"),
                args.getInt("min-mapq", 10),
                args.flag("keep-duplicates"),
                args.get("library", "unstranded"),
                args.getInt("terminal-anchor-tolerance", 5),
                args.getInt("pas-merge-distance", 24),
                args.getInt("ipa-window", 100),
                !args.flag("no-ipa"),
                args.getInt("min-segment-length", 30),
                args.getInt("min-sample-evidence", 10),
                args.getInt("min-informative-replicates", 3),
                args.flag("write-per-sample"),
                args.flag("write-diagnostics"),
                args.getDouble("min-delta-pau", 0.10),
                args.getDouble("min-ale-switch", 0.10),
                args.getDouble("min-ipa-delta", 0.10),
                args.getDouble("fdr", 0.05),
                args.getDouble("wpi-threshold", 0.10)
        );
    }

    private static void runRbpEnrichment(Args args) throws Exception {
        RbpEnrichment.run(
                args.require("sites"),
                args.require("genes"),
                args.require("fasta"),
                args.require("motifs"),
                args.require("out"),
                args.getInt("upstream", 100),
                args.getInt("downstream", 50),
                args.getInt("splice-upstream", 200),
                args.getInt("splice-downstream", 100),
                args.getDouble("site-fdr", 0.05),
                args.getDouble("min-delta-pau", 0.10),
                args.getDouble("motif-score-frac", 0.80),
                args.getInt("background-ratio", 3),
                args.getInt("permutations", 1000),
                Long.parseLong(args.get("seed", "20260830")),
                args.flag("include-no-directional-shift")
        );
    }

    private static void printCommandHelp(String command) {
        switch (command) {
            case "run-direct":
                printDirectHelp(false);
                return;
            case "run-isoquant":
                printDirectHelp(true);
                return;
            case "short-read":
            case "run-short":
                printShortReadHelp();
                return;

            case "sc":
                printScHelp();
                return;

            case "plot-region":
                printPlotRegionHelp();
                return;
            case "rbp-enrich":
                printRbpHelp();
                return;
            case "ends":
                printEndsHelp();
                return;
            case "citation":
            case "--citation":
                Citation.printCitation();
                return;
            case "version":
            case "--version":
            case "-v":
                System.out.println("Usage: java -jar APAURA.jar version");
                return;
            default:
                throw new IllegalArgumentException("Unknown command: " + command);
        }
    }

    private static void printHelp() {
        printLines(
                "Commands:",
                "  ends          Extract single-base transcript-oriented 3' aligned ends",
                "  run-direct    Direct RNA / true 3'-end differential APA + architecture",
                "  run-isoquant  IsoQuant corrected-read differential APA",
                "  short-read    Conventional RNA-seq coverage-based TUTR/ALE/IPA analysis",
                "  run-short     Alias for short-read",
                "  sc            Single-cell APA analysis powered by sc3UTR",
                "  gtf-sites      Build candidate poly(A) sites from a GTF/GFF annotation",
                "  import-sites   Import PAS calls from an external APA caller",
                "  merge-sites    Merge multiple PAS catalogs into a consensus catalog",
                "  quantify       Quantify UMI-deduplicated PAS usage",
                "  reference      Build grouped PAS usage reference tables",
                "  compare        Compare PAS usage between biological conditions",
                "  quantify       Quantify UMI-deduplicated PAS usage from single-cell alignments",
                "  rbp-enrich    RBP motif enrichment around regulated PAS/ALE acceptors",
                "  citation      Print citation",
                "  version       Print version",
                "  plot-region   locus plot for APAURA results",

                "",
                "Long-read/direct required:",
                "  --samples samples.tsv     columns: sample, condition, bam/file/path",
                "  --gtf annotation.gtf[.gz]",
                "  --out DIR",
                "Optional:",
                "  --control control --case treated",
                "  --cluster-distance 20",
                "  --cluster-representative mode|median   [mode]",
                "  --min-gene-reads 10       pooled per condition",
                "  --min-site-reads 5        pooled across selected conditions",
                "  --min-pau 0.01            pooled fraction (0.01 = 1%)",
                "  --min-mapq 10",
                "  --min-delta-pau 0.10 --fdr 0.05 --wpi-threshold 0.10",
                "  --min-internal-exon-delta 0.10 --min-ale-switch 0.10",
                "  --min-ipa-delta 0.10 --min-tutr-shift 0.10",
                "",
                "IsoQuant sample sheet:",
                "  sample\tcondition\tisoquant_dir",
                "  C1\tcontrol\t/path/C1_isoquant",
                "  T1\ttreated\t/path/T1_isoquant",
                "",
                "Conventional short-read required:",
                "  --samples samples.tsv --gtf annotation.gtf --out DIR",
                "Optional:",
                "  --mode annotated",
                "  --library unstranded|forward|reverse",
                "  --terminal-anchor-tolerance 5 --pas-merge-distance 24",
                "  --min-segment-length 30 --min-sample-evidence 10",
                "  --min-mapq 10 --min-delta-pau 0.10 --fdr 0.05 --wpi-threshold 0.10",
                "",
                "RBP enrichment required:",
                "  --sites apa_sites.tsv --genes apa.tsv",
                "  --fasta genome.fa --motifs rbp_motifs.meme --out DIR",
                "Optional:",
                "  --upstream 100 --downstream 50",
                "  --splice-upstream 200 --splice-downstream 100",
                "  --site-fdr 0.05 --min-delta-pau 0.10",
                "  --motif-score-frac 0.80 --background-ratio 3",
                "  --permutations 1000 --seed 20260830",
                "  --include-no-directional-shift",
                "",

                "Examples:",
                "  java -jar APAURA.jar run-direct --samples samples.tsv --gtf gencode.gtf "
                        + "--cluster-distance 0 --out direct_results",
                "  java -jar APAURA.jar run-isoquant --samples isoquant_samples.tsv "
                        + "--gtf gencode.gtf --out isoquant_results",
                "  java -jar APAURA.jar short-read --samples samples.tsv --gtf gencode.gtf "
                        + "--library unstranded --out short_results",
                "  java -jar APAURA.jar rbp-enrich --sites results/apa_sites.tsv "
                        + "--genes results/apa.tsv --fasta genome.fa --motifs rbp.meme "
                        + "--out rbp_results",
                "",
                "For command-specific help:",
                "  java -jar APAURA.jar <command> --help"
        );
    }

    private static void printScHelp() {
        printLines(
                "Usage:",
                "  java -jar APAURA.jar sc <command> [options]",
                "",
                "Single-cell alternative polyadenylation analysis",
                "",
                "Commands:",
                "  gtf-sites      Build candidate poly(A) sites from a GTF/GFF annotation",
                "  import-sites   Import PAS calls from an external APA caller",
                "  merge-sites    Merge multiple PAS catalogs into a consensus catalog",
                "  quantify       Quantify UMI-deduplicated PAS usage",
                "  reference      Build grouped PAS usage reference tables",
                "  compare        Compare PAS usage between biological conditions",
                "",
                "Run:",
                "  java -jar APAURA.jar sc <command> --help"
        );
    }

    private static void printScImportSitesHelp() {
        printLines(
                "Usage:",
                "  java -jar APAURA.jar sc import-sites --input caller.tsv --source NAME --out sites.tsv [options]",
                "",
                "Import candidate poly(A) sites from an external APA caller.",
                "",
                "Required:",
                "  --input FILE            input caller table",
                "  --source NAME           source name, for example SCAPE",
                "  --out FILE              output PAS catalog",
                "",
                "Options:",
                "  --format FORMAT         auto|scape|sierra|scapatrap|bed [auto]",
                "  --*-column VALUE        override detected input columns",
                "",
                "Example:",
                "  java -jar APAURA.jar sc import-sites --input sites.bed --source custom --format bed --out sites.tsv"
        );
    }

    private static void printScQuantifyHelp() {
        printLines(
                "Usage:",
                "  java -jar APAURA.jar sc quantify --sites sites.tsv --alignment sample.sam|bam|cram",
                "    --sample-id ID --out counts.tsv [options]",
                "",
                "Quantify single-cell poly(A) site usage from SAM, BAM, or CRAM alignments.",
                "",
                "Required:",
                "  --sites FILE            PAS reference",
                "  --alignment FILE        SAM, BAM, or CRAM alignment",
                "  --sample-id ID          biological sample identifier",
                "  --out FILE              output PAS count table",
                "",
                "Options:",
                "  --barcode-tag TAG       cell barcode tag [CB]",
                "  --umi-tag TAG           UMI tag [UB]",
                "  --gene-tag TAG          gene assignment tag [GX]",
                "  --gene-tag-mode MODE    plain|bd-xf [plain]",
                "  --site-window INT       maximum distance to compatible PAS [50]",
                "  --three-prime-by MODE   gene|read|opposite-read|nearest-end [gene]",
                "",
                "Examples:",
                "  java -jar APAURA.jar sc quantify --sites sites.tsv --alignment sample.bam \\",
                "    --sample-id S1 --out S1.counts.tsv",
                "",
                "  java -jar APAURA.jar sc quantify --sites sites.tsv --alignment sample.bam \\",
                "    --sample-id S1 --gene-tag XF --gene-tag-mode bd-xf --out S1.counts.tsv"
        );
    }

    private static void printScReferenceHelp() {
        printLines(
                "Usage:",
                "  java -jar APAURA.jar sc reference --sites sites.tsv --counts counts.tsv [options] --out DIR",
                "",
                "Build grouped single-cell poly(A) site usage reference tables.",
                "",
                "Required:",
                "  --sites FILE            PAS reference",
                "  --counts FILE           quantified PAS counts",
                "  --out DIR               output directory",
                "",
                "Options:",
                "  --cells FILE            cell metadata table",
                "  --group-by COLUMNS      comma-separated grouping columns",
                "",
                "Example:",
                "  java -jar APAURA.jar sc reference --sites sites.tsv --counts counts.tsv \\",
                "    --cells cells.tsv --group-by cell_type --out reference_results"
        );
    }

    private static void printScCompareHelp() {
        printLines(
                "Usage:",
                "  java -jar APAURA.jar sc compare --sites sites.tsv --samples samples.tsv \\",
                "    --counts counts.tsv --control NAME --treated NAME --out DIR [options]",
                "",
                "Compare single-cell poly(A) site usage between biological conditions.",
                "",
                "Required:",
                "  --sites FILE            PAS reference",
                "  --samples FILE          biological sample metadata",
                "  --counts FILE           quantified PAS counts",
                "  --control NAME          control condition",
                "  --treated NAME          treated condition",
                "  --out DIR               output directory",
                "",
                "Cell grouping:",
                "  --cells FILE            cell metadata table",
                "  --group-by COLUMNS      e.g. cell_type,cell_subtype",
                "  --group-column COLUMN   legacy single-column grouping",
                "  --group-separator TEXT  separator for combined group labels [::]",
                "  --label-group-columns   label combined group columns",
                "",
                "Thresholds:",
                "  --min-sample-count INT      [5]",
                "  --min-condition-count INT   [20]",
                "  --min-delta FLOAT            [0.10]",
                "  --min-jsd FLOAT              [0.10]",
                "  --fdr FLOAT                  [0.05]",
                "  --distribution-permutations INT [999]",
                "  --min-interaction FLOAT      [0.10]",
                "  --no-interactions",
                "",
                "Interaction = treated(group1-group2) - control(group1-group2)"
        );
    }

    private static void printScMergeSitesHelp() {
        printLines(
                "Usage:",
                "  java -jar APAURA.jar sc merge-sites --inputs a.tsv,b.tsv --out consensus.tsv [options]",
                "",
                "Merge multiple poly(A) site catalogs into a consensus PAS catalog.",
                "",
                "Required:",
                "  --inputs FILES          comma-separated PAS catalogs",
                "  --out FILE              output consensus PAS catalog",
                "",
                "Options:",
                "  --sources NAMES         comma-separated source names, for example GTF,SCAPE",
                "  --cluster-bp INT        cluster nearby PAS positions [24]",
                "  --min-sources INT       minimum number of supporting sources [1]",
                "",
                "Example:",
                "  java -jar APAURA.jar sc merge-sites --inputs gtf.tsv,scape.tsv \\",
                "    --sources GTF,SCAPE --cluster-bp 24 --out consensus.tsv"
        );
    }

    private static void printScGtfSitesHelp() {
        printLines(
                "Usage:",
                "  java -jar APAURA.jar sc gtf-sites --annotation annotation.gtf|gff3 --out sites.tsv [options]",
                "",
                "Build a strand-aware candidate 3' end / poly(A) site catalog",
                "from transcript annotations.",
                "",
                "Required:",
                "  --annotation FILE       GTF or GFF3 transcript annotation",
                "  --out FILE              output PAS catalog",
                "",
                "Options:",
                "  --cluster-bp INT        cluster nearby annotated transcript ends [24]",
                "",
                "Example:",
                "  java -jar APAURA.jar sc gtf-sites --annotation gencode.gtf --out sites.tsv"
        );
    }
    private static void printPlotRegionHelp() {

        printLines(
                "Usage:",
                "  java -jar APAURA.jar plot-region [--results DIR] (--gene NAME[,NAME...] | --genes-file FILE | --region chr:start-end) [options]",
                "",
                "Target:",
                "  --gene NAME[,NAME...]      gene symbol/Ensembl ID; comma-separated genes allowed",
                "  --genes-file FILE          gene symbols/Ensembl IDs; comments with # allowed",
                "  --region REGION            chr:start-end; exclusive with gene targets",
                "",
                "Input:",
                "  --results DIR              optional completed APAURA result directory",
                "                             with results: coverage + tested TUTR/ALE/IPA",
                "                             without results: standalone short-read locus plot",
                "  --samples FILE             sample sheet; required standalone",
                "  --gtf FILE                 annotation GTF/GTF.GZ; required standalone",
                "  --control NAME             control condition; required standalone",
                "  --case NAME                case condition; required standalone",
                "  --mode MODE                short-read|run-direct|run-isoquant",
                "                             standalone defaults to short-read",
                "  --library MODE             unstranded|forward|reverse [unstranded]",
                "",
                "Plot options:",
                "  --flank INT                gene flank [1000]",
                "  --merge mean|median        replicate merge for short-read coverage [mean]",
                "  --show-replicates          overlay individual replicate tracks",
                "  --show-overview            add whole-gene coverage above the 3′ terminal zoom",
                "  --verbose                  show per-sample coverage/cache diagnostics",
                "  --bin-size INT             short-read coverage bin size [10]",
                "  --min-mapq INT             minimum mapping quality [10]",
                "  --samtools FILE            samtools executable [samtools]",
                "  --keep-duplicates          retain duplicate reads",
                "  --no-ipa                   omit IPA panel when results are available",
                "  --genomic-orientation      preserve genomic coverage orientation",
                "  --format pdf|png|svg|all   output format [all]",
                "  --width FLOAT              figure width [12]",
                "  --height FLOAT             figure height [9]",
                "  --rscript FILE             Rscript executable [Rscript]",
                "  --out PREFIX|DIR           prefix for one target; directory for multiple genes",
                "  --show-transcripts         show individual annotated transcript models",
                "  --max-transcripts INT      maximum transcript/terminal models displayed [20]",
                "",
                "Standalone behavior:",
                "  A completed APAURA result directory is not required.",
                "  Coverage and annotation are drawn directly from BAM/SAM + GTF.",
                "  Differential TUTR/IPA/ALE panels are absent because no APA test was run.",
                "",
                "Examples:",
                "  # Completed APAURA results",
                "  java -jar APAURA.jar plot-region --results results --gene HNRNPA2B1",
                "",
                "  # Standalone gene",
                "  java -jar APAURA.jar plot-region --samples samples.tsv --gtf gencode.gtf \\",
                "    --gene HNRNPA2B1 --control WT --case MUT --library reverse --show-overview",
                "",
                "  # Standalone genomic region",
                "  java -jar APAURA.jar plot-region --samples samples.tsv --gtf gencode.gtf \\",
                "    --region chr7:26170000-26197000 --control WT --case MUT --library reverse"
        );
    }

    private static void printDirectHelp(boolean isoquant) {
        String command = isoquant ? "run-isoquant" : "run-direct";
        String sampleColumn = isoquant ? "isoquant_dir" : "bam/file/path";

        printLines(
                "Usage:",
                "  java -jar APAURA.jar " + command
                        + " --samples samples.tsv --gtf annotation.gtf --out DIR [options]",
                "",
                "Required:",
                "  --samples FILE            sample sheet: sample, condition, " + sampleColumn,
                "  --gtf FILE                GTF/GTF.GZ annotation",
                "  --out DIR                 output directory",
                "",
                "Options:",
                "  --control NAME            control condition [control]",
                "  --case NAME               case condition [treated]",
                "  --cluster-distance INT    PAS clustering distance [20]",
                "  --cluster-representative mode|median [mode]",
                "  --min-gene-reads INT      pooled reads per condition [10]",
                "  --min-site-reads INT      pooled PAS reads [5]",
                "  --min-pau FLOAT           minimum pooled PAS usage [0.01]",
                "  --min-mapq INT            minimum mapping quality [10]",
                "  --keep-duplicates         retain duplicate-marked reads",
                "  --strand-mode MODE        alignment strand mode [alignment]",
                "  --annotation-slop INT     downstream annotation window [1000]",
                "  --min-delta-pau FLOAT     PAS/branch usage effect threshold [0.10]",
                "  --min-ale-switch FLOAT    minimum ALE branch-distribution switch [0.10]",
                "  --fdr FLOAT               FDR threshold [0.05]",
                "  --wpi-threshold FLOAT     shortening/lengthening threshold [0.10]",
                "  --min-internal-exon-delta FLOAT [0.10]",
                "  --min-ale-switch FLOAT    ALE switch threshold [0.10]",
                "  --min-ipa-delta FLOAT     IPA gain/loss threshold [0.10]",
                "  --min-tutr-shift FLOAT    TUTR shift threshold [0.10]",
                "",
                "Example:",
                "  java -jar APAURA.jar " + command
                        + " --samples samples.tsv --gtf gencode.gtf --out results"
        );
    }

    private static void printShortReadHelp() {
        printLines(
                "Usage:",
                "  java -jar APAURA.jar short-read --samples samples.tsv --gtf annotation.gtf --out DIR [options]",
                "",
                "Required:",
                "  --samples FILE            sample sheet: sample, condition, bam/file/path",
                "  --gtf FILE                GTF/GTF.GZ annotation",
                "  --out DIR                 output directory",
                "",
                "Options:",
                "  --mode annotated          annotation-guided TUTR/ALE/IPA mode [annotated]",
                "  --control NAME            control condition [control]",
                "  --case NAME               case condition [treated]",
                "  --library MODE            unstranded|forward|reverse [unstranded]",
                "  --min-mapq INT            minimum mapping quality [10]",
                "  --keep-duplicates         retain duplicate-marked reads",
                "  --terminal-anchor-tolerance INT [5]",
                "  --pas-merge-distance INT  merge nearby annotated PAS [24]",
                "  --min-segment-length INT  minimum coverage segment length [30]",
                "  --min-sample-evidence INT minimum sample evidence [10]",
                "  --min-informative-replicates INT  informative biological replicates required per condition [3]",
                "  --write-per-sample        write per-sample usage tables under details/",
                "  --write-diagnostics       write developer/full-detail tables under details/",
                "  --min-delta-pau FLOAT     PAS/branch usage effect threshold [0.10]",
                "  --min-ale-switch FLOAT    minimum ALE branch-distribution switch [0.10]",
                "  --fdr FLOAT               FDR threshold [0.05]",
                "  --ipa-window INT          IPA quantification window [100]",
                "  --min-ipa-delta FLOAT     IPA gain/loss threshold [0.10]",
                "  --wpi-threshold FLOAT     shortening/lengthening threshold [0.10]",
                "",
                "Primary outputs:",
                "  apa.tsv                   integrated one-row-per-gene APA result",
                "  apa_events.tsv            TUTR/ALE/IPA event-level differential results",
                "  apa_sites.tsv             terminal PAS-level differential results",
                "  details/                  architecture and mechanism-specific supporting tables",
                "",
                "Example:",
                "  java -jar APAURA.jar short-read --samples samples.tsv --gtf gencode.gtf "
                        + "--library unstranded --out short_results"
        );
    }

    private static void printRbpHelp() {
        printLines(
                "Usage:",
                "  java -jar APAURA.jar rbp-enrich --sites apa_sites.tsv --genes apa.tsv "
                        + "--fasta genome.fa --motifs rbp.meme --out DIR [options]",
                "",
                "Required:",
                "  --sites FILE              APAURA apa_sites.tsv",
                "  --genes FILE              APAURA apa.tsv",
                "  --fasta FILE              reference genome FASTA",
                "  --motifs FILE             RBP motifs in MEME format",
                "  --out DIR                 output directory",
                "",
                "Options:",
                "  --upstream INT            PAS upstream window [100]",
                "  --downstream INT          PAS downstream window [50]",
                "  --splice-upstream INT     ALE acceptor upstream window [200]",
                "  --splice-downstream INT   ALE acceptor downstream window [100]",
                "  --site-fdr FLOAT          APA-site FDR filter [0.05]",
                "  --min-delta-pau FLOAT     minimum absolute delta PAU [0.10]",
                "  --motif-score-frac FLOAT  normalized motif threshold [0.80]",
                "  --background-ratio INT    matched background sites per foreground [3]",
                "  --permutations INT        permutation count [1000]",
                "  --seed INT                random seed [20260830]",
                "  --include-no-directional-shift",
                "",
                "Example:",
                "  java -jar APAURA.jar rbp-enrich --sites results/apa_sites.tsv "
                        + "--genes results/apa.tsv --fasta genome.fa --motifs rbp.meme "
                        + "--out rbp_results"
        );
    }

    private static void printEndsHelp() {
        printLines(
                "Usage:",
                "  java -jar APAURA.jar ends --bam sample.bam --out ends.bed [options]",
                "",
                "Required:",
                "  --bam FILE                BAM/SAM alignment file",
                "  --out FILE                output single-base BED-like table",
                "",
                "Options:",
                "  --min-mapq INT            minimum mapping quality [10]",
                "  --keep-duplicates         retain duplicate-marked reads",
                "  --strand-mode MODE        alignment strand mode [alignment]"
        );
    }

    private static void printLines(String... lines) {
        for (String line : lines) {
            System.out.println(line);
        }
    }
}
