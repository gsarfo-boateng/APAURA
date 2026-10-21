# Development guide

## Full-build rule

After every source-code change:

```bash
mvn clean package
./build.sh
```

Do not rely only on compiling individual Java files.

Then smoke-test:

```bash
apaura --version
apaura short-read --help
apaura plot-region --help
```

## R resource validation

After editing the bundled region renderer:

```bash
Rscript -e 'parse(file="src/main/resources/org/apaura/plot/APAURA_region_plot.R"); cat("R syntax OK\n")'
```

Then rebuild the complete JAR so the updated resource is embedded.

## Tests

```bash
mvn test
mvn clean package
./tests/run_tests.sh
```

## Packaging checks

```bash
jar tf APAURA.jar | grep 'org/apaura/stats/MultipleTesting.class'
jar tf APAURA.jar | grep 'org/apaura/analysis/ShortReadDifferential.class'
```

## Regression expectations

Short-read tests should cover:

- adequate replicate support;
- sparse replicate support;
- genuine 1-vs-1 Fisher behavior;
- BH correction only over finite tested p-values;
- observed vs inferential direction;
- IPA replicate support;
- compact output gating.

Plot tests should cover:

- plus and minus strands;
- fixed genomic x coordinates;
- strand-aware extraction;
- display-only mirroring;
- no replicate-color leakage between facets;
- normalization cache miss then hit;
- multi-gene output;
- `--genes-file` handling.
