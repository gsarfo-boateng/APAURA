# APAURA `plot-region` patch

This patch adds a unified Java command that delegates publication-quality locus
rendering to a bundled R script.

## New command

```bash
java -jar APAURA.jar plot-region \
  --results RESULTS_DIR \
  --gene HNRNPA2B1
```

The command supports `short-read`, `run-direct`, and `run-isoquant`.

- `short-read`: merged normalized BAM coverage by condition + ΔPAU + optional ΔIPA.
- `run-direct`: direct-RNA PAS/PAU evidence + ΔPAU + optional IPA.
- `run-isoquant`: corrected-read PAS/PAU evidence + ΔPAU + optional IPA.
- Negative-strand genes are automatically displayed in transcript 5'→3'
  orientation. Use `--genomic-orientation` to retain increasing genomic
  coordinates left-to-right.

## 1. Copy the three new files

From this patch into the project:

```text
src/main/java/org/apaura/plot/PlotRegionCommand.java
src/main/java/org/apaura/util/RunMetadata.java
src/main/resources/org/apaura/plot/APAURA_region_plot.R
```

Maven automatically packages `src/main/resources`.

## 2. Edit the REAL source `src/main/java/org/apaura/APAURA.java`

Do not edit a decompiled `target/classes` or Fernflower copy.

Add imports:

```java
import java.io.File;
import java.util.LinkedHashMap;
import org.apaura.plot.PlotRegionCommand;
import org.apaura.util.RunMetadata;
```

### Add the command to `runCommand`

```java
case "plot-region":
    PlotRegionCommand.run(args);
    return;
```

Use your real local variable name in place of `args`.

### Add command help dispatch

```java
case "plot-region":
    printPlotRegionHelp();
    return;
```

### Add to global command list

```text
plot-region   Publication-style locus plot for short-read/direct/IsoQuant results
```

### Add this help method

```java
private static void printPlotRegionHelp() {
    printLines(
        "Usage:",
        "  java -jar APAURA.jar plot-region --results DIR (--gene NAME | --region chr:start-end) [options]",
        "",
        "Required:",
        "  --results DIR             completed APAURA result directory",
        "  --gene NAME               gene symbol/Ensembl ID (exclusive with --region)",
        "  --region REGION           chr:start-end (exclusive with --gene)",
        "",
        "Normally read from apaura_run.tsv:",
        "  --samples FILE            override sample sheet",
        "  --gtf FILE                override GTF",
        "  --control NAME            override control condition",
        "  --case NAME               override case condition",
        "  --mode MODE               short-read|run-direct|run-isoquant",
        "",
        "Plot options:",
        "  --flank INT               gene flank [1000]",
        "  --merge mean|median       replicate merge for short-read coverage [mean]",
        "  --show-replicates         overlay individual short-read replicate tracks",
        "  --bin-size INT            short-read coverage bin size [10]",
        "  --min-mapq INT            short-read minimum MAPQ [run metadata or 10]",
        "  --samtools FILE           samtools executable [samtools]",
        "  --no-ipa                  hide IPA panel",
        "  --genomic-orientation     do not flip negative-strand genes",
        "  --min-delta-pau FLOAT     APA effect guide [0.10]",
        "  --min-ipa-delta FLOAT     IPA effect guide [0.10]",
        "  --fdr FLOAT               significance threshold [0.05]",
        "  --format pdf|png|svg|all  output format [all]",
        "  --width FLOAT             figure width [12]",
        "  --height FLOAT            figure height [9]",
        "  --rscript FILE            Rscript executable [Rscript]",
        "  --out PREFIX              output prefix [DIR/plots/TARGET.APAURA]",
        "",
        "Example:",
        "  java -jar APAURA.jar plot-region --results short_results --gene HNRNPA2B1"
    );
}
```

## 3. Write run metadata after successful analyses

The metadata belongs in the CLI wrapper because the wrapper already knows the
sample sheet, GTF, conditions, and mode. This avoids changing analysis-engine
method signatures.

### `run-direct`

After `DirectApaPipeline.run(...)` succeeds:

```java
LinkedHashMap<String, String> meta = new LinkedHashMap<String, String>();
meta.put("input_format", inputFormat);
meta.put("min_mapq", String.valueOf(minMapq));
meta.put("keep_duplicates", String.valueOf(keepDuplicates));

RunMetadata.write(
    new File(outDir),
    VERSION,
    "run-direct",
    sampleSheet,
    gtf,
    control,
    caseCondition,
    meta
);
```

For the IsoQuant branch, write analysis as `"run-isoquant"` rather than
`"run-direct"`.

A compact pattern is:

```java
String analysis = inputFormat.equalsIgnoreCase("isoquant")
        ? "run-isoquant"
        : "run-direct";
```

### `short-read`

After `ShortReadPipeline.run(...)` succeeds:

```java
LinkedHashMap<String, String> meta = new LinkedHashMap<String, String>();
meta.put("library", library);
meta.put("min_mapq", String.valueOf(minMapq));
meta.put("keep_duplicates", String.valueOf(keepDuplicates));
meta.put("pas_merge_distance", String.valueOf(pasMergeDistance));
meta.put("ipa_enabled", String.valueOf(runIpa));

RunMetadata.write(
    new File(outDir),
    VERSION,
    "short-read",
    sampleSheet,
    gtf,
    control,
    caseCondition,
    meta
);
```

Important: store the CLI values in local variables before calling the pipeline,
then reuse those variables in `RunMetadata.write`. Do not repeatedly parse flags.

## 4. Legacy result directories

Existing directories do not yet contain `apaura_run.tsv`. Plot them by supplying
the missing values explicitly:

```bash
java -jar APAURA.jar plot-region \
  --results /path/old_results \
  --mode short-read \
  --samples /path/samples.tsv \
  --gtf /path/annotation.gtf \
  --control U1_WT \
  --case U1_MUT \
  --samtools /home/george/miniconda3/bin/samtools \
  --gene HNRNPA2B1
```

## 5. R dependencies

The bundled renderer needs:

```r
ggplot2
dplyr
readr
stringr
purrr
tibble
patchwork
rtracklayer
GenomicRanges
```

The renderer does NOT require `optparse`; Java owns CLI parsing.

## 6. Full validation

APAURA development rule: run the entire build after every code change.

```bash
cd /media/george/Extreme_SSD/APAURA

mvn clean package
./build.sh

java -jar APAURA.jar --version
java -jar APAURA.jar plot-region --help
```

Then validate the current short-read result:

```bash
java -jar APAURA.jar plot-region \
  --results /media/george/George/CLL_Map_fastq/apaura_short_results_IPA \
  --mode short-read \
  --samples /media/george/George/CLL_Map_fastq/apaura_short_results_IPA/samples.tsv \
  --gtf /media/george/George/Erica/Human_Lymphiod/qapa/reference/gencode.v49.primary_assembly.annotation.gtf \
  --control U1_WT \
  --case U1_MUT \
  --samtools /home/george/miniconda3/bin/samtools \
  --gene HNRNPA2B1 \
  --out /media/george/George/CLL_Map_fastq/apaura_short_results_IPA/HNRNPA2B1_APAURA
```

Expected behavior for HNRNPA2B1:
- U1_WT and U1_MUT merged coverage tracks.
- negative strand automatically flipped to transcript 5'→3' orientation.
- ΔPAU terminal PAS panel.
- ΔIPA panel when `ipa_sites.tsv` exists.
- gene model arrow follows transcription.
