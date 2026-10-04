# APAURA figure generation

APAURA includes `scripts/plot_apaura.R`, a companion R script that converts the standard **PAS → PAU → APA → RBP** result files into publication-ready summary figures.

The plotting script is deliberately separate from the Java analysis engine. APAURA itself does not require R; R is only needed when figures are requested.

## Requirements

R 4.1 or newer is recommended. The script uses:

- `ggplot2`
- `dplyr`
- `tidyr`
- `readr`
- `scales`
- `tibble`

Install once with:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "scales", "tibble"))
```

## Basic usage

After a APAURA run:

```bash
Rscript scripts/plot_apaura.R \
  --input apaura_results
```

By default figures are written to:

```text
apaura_results/figures/
```

If RBP enrichment was written to a separate directory:

```bash
Rscript scripts/plot_apaura.R \
  --input apaura_results \
  --rbp rbp_results \
  --out apaura_results/figures
```

## Options

```text
--input DIR          APAURA results directory
--rbp DIR            RBP results directory
--out DIR            output figure directory
--fdr FLOAT          FDR threshold; default 0.05
--delta-pau FLOAT    absolute delta-PAU threshold; default 0.10
--top INT            number of top genes/RBPs shown; default 20
--width FLOAT        default PDF/PNG width; default 8 inches
--height FLOAT       default PDF/PNG height; default 6 inches
```

Example:

```bash
Rscript scripts/plot_apaura.R \
  --input results \
  --rbp results/rbp \
  --fdr 0.05 \
  --delta-pau 0.10 \
  --top 25 \
  --out results/figures
```

## Output

Each applicable plot is written as both PDF and 300-dpi PNG. The script also creates:

```text
APAURA_figures.pdf
figure_manifest.tsv
plot_run_info.txt
```

`APAURA_figures.pdf` is a multi-page figure book containing all plots that could be generated from the available outputs.

The script is schema-aware. A short-read analysis, for example, does not contain all long-read IPA/ALE columns; plots requiring unavailable columns are skipped rather than failing.

## Figure families

Depending on the available results, the script can generate:

### PAS

- PAS counts by genomic region
- PAS counts by architecture role
- number of PAS per gene
- PAS-supporting reads per sample
- short-read observed-versus-fitted coverage

### PAU

- PAU distribution by condition
- sample-to-sample PAU correlation heatmap
- PCA of sample PAU profiles

### APA site-level

- delta-PAU/FDR volcano plot
- control versus case PAU scatter plot
- delta PAU by PAS architecture role

### APA gene-level

- APA architecture counts
- APA direction counts
- architecture-by-direction composition
- control versus case WPI
- delta WPI by architecture
- TUTR WPI shift
- ALE WPI shift
- IPA usage shift
- termination-shift index versus Jensen-Shannon divergence
- mechanism-specific TUTR/ALE/IPA call matrix
- top architecture-aware APA effect sizes

### RBP

- top RBP motif-enrichment dot plot
- foreground versus matched-background motif prevalence
- positional motif profiles relative to PAS/ALE splice acceptors
- RBP-by-APA-mechanism enrichment heatmap

## Interpretation

The plotting script does not recalculate the primary APAURA statistics. It visualizes values already present in `pau.tsv`, `apa_sites.tsv`, `apa.tsv`, and `rbp.tsv`.

For TUTR events, negative positional shifts indicate shortening and positive shifts indicate lengthening. For ALE, branch shifts should be interpreted using `ale_direction`. IPA should be interpreted using `ipa_direction` and condition-level IPA usage rather than as terminal-exon shortening or lengthening.

## Reproducibility

`plot_run_info.txt` records plotting thresholds, paths, and the R version. `figure_manifest.tsv` lists every generated figure and a short description, making it easier to include selected figures in manuscripts or reports.
