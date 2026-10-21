# Release checklist

## Code

- [ ] Intended source changes committed.
- [ ] No debug code remains.
- [ ] No private workstation paths are hard-coded.
- [ ] Version and citation text agree.

## Full build

```bash
mvn clean package
./build.sh
```

- [ ] Zero unit-test failures.
- [ ] `./build.sh` succeeds.

## JAR contents

```bash
jar tf APAURA.jar | grep 'org/apaura/stats/MultipleTesting.class'
jar tf APAURA.jar | grep 'org/apaura/analysis/ShortReadDifferential.class'
```

- [ ] Required classes present.

## CLI

```bash
apaura --version
apaura --help
apaura short-read --help
apaura run-direct --help
apaura run-isoquant --help
apaura plot-region --help
apaura rbp-enrich --help
```

- [ ] Help matches actual behavior.
- [ ] Multi-gene plot syntax documented.
- [ ] `--library` documented correctly.
- [ ] Plot help states genomic x coordinates are preserved.

## R renderer

```bash
Rscript -e 'parse(file="src/main/resources/org/apaura/plot/APAURA_region_plot.R"); cat("R syntax OK\n")'
```

- [ ] Plus-strand plot inspected.
- [ ] Minus-strand plot inspected.
- [ ] Replicate colors do not mix.
- [ ] Normalization cache miss/hit verified.
- [ ] Multi-gene plotting verified.

## Integration tests

```bash
./tests/run_tests.sh
```

- [ ] Short-read smoke test.
- [ ] Direct-mode smoke test.
- [ ] IsoQuant smoke test.
- [ ] RBP-enrichment smoke test.

## Documentation

- [ ] README matches current CLI.
- [ ] Output schema matches writers.
- [ ] CHANGELOG updated.
- [ ] LICENSE present.
- [ ] CITATION.cff present.

## GitHub release

- [ ] Tag release.
- [ ] Attach final `APAURA.jar`.
- [ ] Attach checksums if distributed.
- [ ] Publish release notes.
