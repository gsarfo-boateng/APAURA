# Validation

## Build

```bash
mvn clean package
./build.sh
```

Require zero failures.

## JAR packaging

```bash
jar tf APAURA.jar | grep 'org/apaura/stats/MultipleTesting.class'
jar tf APAURA.jar | grep 'org/apaura/analysis/ShortReadDifferential.class'
```

## CLI smoke tests

```bash
java -jar APAURA.jar --version
java -jar APAURA.jar --help
java -jar APAURA.jar short-read --help
java -jar APAURA.jar run-direct --help
java -jar APAURA.jar run-isoquant --help
java -jar APAURA.jar plot-region --help
java -jar APAURA.jar rbp-enrich --help
```

## Integration tests

```bash
./tests/run_tests.sh
```

## Replicate-support regression

A replicated event with large pooled effect but too few informative replicates should return:

```text
status = insufficient_replicate_support
pvalue = NA
fdr = NA
```

It must not be rescued by pooled Fisher testing.

## Plot regression

Inspect at least one plus-strand and one minus-strand gene. Confirm:

- condition colors remain separate;
- genomic coordinates align across all panels;
- minus-strand mirroring is display-only;
- PAS/IPA coordinates match result tables;
- cache miss occurs only when normalization must be computed;
- later plots use cache hits;
- multi-gene plotting writes one plot per target.
