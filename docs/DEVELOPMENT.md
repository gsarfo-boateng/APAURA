# APAURA development notes

## Project structure

APAURA uses the standard Maven layout:

```text
src/main/java/org/apaura/      production Java source
src/test/java/org/apaura/      JUnit unit tests
tests/                         end-to-end synthetic regression tests
examples/                      small example/test datasets
scripts/                       companion analysis and plotting scripts
docs/                          algorithm and user documentation
```

Generated build products belong in `target/` or `build/` and are excluded from version control.

## Development JDK

The project targets Java 8-compatible bytecode. JDK 17 or 21 is suitable for development.

```bash
mvn clean package
```

runs unit tests and builds `target/APAURA.jar`.

## Test layers

APAURA uses two complementary test layers.

### JUnit tests

`src/test/java/org/apaura/` tests small mathematical and parsing units, including:

- Fisher exact testing;
- Benjamini-Hochberg correction;
- weighted isotonic regression;
- tie-aware ECDF/KS statistics;
- beta-binomial differential usage;
- PAS clustering;
- CIGAR parsing.

Run:

```bash
mvn test
```

### Synthetic integration tests

`tests/run_tests.sh` creates synthetic datasets, builds APAURA, executes complete workflows, and checks expected output values and file concordance.

```bash
./tests/run_tests.sh
```

When refactoring algorithms, the integration outputs should remain unchanged unless the scientific method or output schema is intentionally changed.

## BAM/BGZF reader scope

APAURA includes `BamLiteReader` and `BgzfInputStream` to support dependency-light sequential BAM streaming.

The reader supports the information APAURA currently needs:

- reference identifier and alignment position;
- mapping quality;
- SAM flags;
- read/query length;
- CIGAR-derived aligned reference blocks.

It does not aim to replace a complete BAM library. In particular, it does not provide indexed random access, optional-tag parsing, pileup APIs, CRAM decoding, or general BAM editing.

The custom reader should remain covered by synthetic BAM/SAM concordance tests. If APAURA later requires broader alignment-format functionality, replacing or supplementing it with a mature external library such as HTSJDK should be considered.

## Code style

Production Java should favor:

- descriptive variable and method names;
- small private helper methods;
- explicit imports;
- technical Javadocs where the implementation is non-obvious;
- line lengths that remain readable in IntelliJ;
- biological/methodological discussion in `docs/`, not inside implementation comments.

## Release checklist

Before a release:

```bash
mvn clean package
./tests/run_tests.sh
java -jar target/APAURA.jar --help
java -jar target/APAURA.jar --citation
```

Also confirm that generated directories (`target/`, `build/`, `tests/out_*`) and macOS metadata files (`._*`, `.DS_Store`) are not included in the source archive.
