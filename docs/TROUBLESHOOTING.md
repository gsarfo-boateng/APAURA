# Troubleshooting

## `NoClassDefFoundError: org/apaura/stats/MultipleTesting`

This is a packaging/build problem, not a biological-data problem.

```bash
find src/main/java -name 'MultipleTesting.java' -print
find target/classes -path '*org/apaura/stats/MultipleTesting.class' -print
jar tf target/APAURA.jar | grep 'org/apaura/stats/MultipleTesting.class'
jar tf APAURA.jar | grep 'org/apaura/stats/MultipleTesting.class'
```

Then run:

```bash
mvn clean package
./build.sh
```

## `plot-region` is slow on the first gene

The first short-read plot may populate `plot_normalization.tsv` by scanning each BAM once. Later plots should use cache hits and be much faster.

## Control and case colors appear in both facets

Use the current renderer. Replicate-level data must contain the same condition-specific facet label used by the merged coverage data.

## Negative-strand coverage is below zero

This is a display convention. Coverage magnitude remains biological coverage; only the plotting sign is changed.

Use `--genomic-orientation` to keep it upright.

## R parser errors

```bash
Rscript -e 'parse(file="src/main/resources/org/apaura/plot/APAURA_region_plot.R"); cat("R syntax OK\n")'
```

Run this before rebuilding after editing the R resource.

## Large effect but `insufficient_replicate_support`

APAURA intentionally does not substitute pooled sequencing depth for biological replication. Inspect the total/any/informative replicate-count columns.

## `min_fdr = NA`

This can indicate insufficient replicate support, insufficient evidence, or no valid finite site-level test. Do not reinterpret `NA` as FDR=1.

## BAM plotting errors

```bash
samtools quickcheck sample.bam
samtools index sample.bam
```

Also verify chromosome naming agrees between BAM and GTF.
