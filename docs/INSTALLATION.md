# Installation

## Core requirements

- Java 8 or newer
- Maven for reproducible source builds
- GTF/GTF.GZ annotation
- coordinate-sorted BAM/SAM where applicable

`plot-region` additionally requires R, `samtools`, indexed BAMs, and the R packages `ggplot2`, `dplyr`, `readr`, `stringr`, `purrr`, `tibble`, `patchwork`, `rtracklayer`, and `GenomicRanges`.

## Build

```bash
git clone <YOUR-APAURA-REPOSITORY-URL>
cd APAURA
mvn clean package
./build.sh
```

APAURA development should always use the full build workflow after source changes.

## Verify the final JAR

```bash
java -jar APAURA.jar --version
java -jar APAURA.jar short-read --help
java -jar APAURA.jar plot-region --help
```

Also verify required runtime classes are packaged:

```bash
jar tf APAURA.jar | grep 'org/apaura/stats/MultipleTesting.class'
jar tf APAURA.jar | grep 'org/apaura/analysis/ShortReadDifferential.class'
```

This catches cases where code compiled against a class that is absent from the final JAR.

## Plotting dependencies

```r
install.packages(c("ggplot2","dplyr","readr","stringr","purrr","tibble","patchwork"))
if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("rtracklayer","GenomicRanges"))
```

Check `samtools`:

```bash
samtools --version
```

Index BAMs when needed:

```bash
samtools index sample.bam
```

## Full validation

```bash
mvn clean package
./build.sh
./tests/run_tests.sh
java -jar APAURA.jar --version
```
