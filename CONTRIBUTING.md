# Contributing to APAURA

Contributions are welcome.

Before a pull request:

1. add/update tests;
2. update CLI/schema documentation when behavior changes;
3. run the complete build;
4. run integration tests.

```bash
mvn clean package
./build.sh
./tests/run_tests.sh
```

For plot changes:

```bash
Rscript -e 'parse(file="src/main/resources/org/apaura/plot/APAURA_region_plot.R")'
```

Development principles:

- preserve Java 8-compatible bytecode;
- keep technology-specific evidence models explicit;
- do not silently substitute sequencing depth for biological replication;
- write unavailable inferential values as `NA`;
- do not treat conventional short-read ends as direct cleavage evidence;
- document primary output-schema changes.
