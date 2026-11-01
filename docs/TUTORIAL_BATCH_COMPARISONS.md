# Tutorial: batch comparisons

For studies with many pairwise contrasts, keep one sample sheet and one result directory per comparison.

```bash
#!/usr/bin/env bash
set -euo pipefail

APAURA="/path/APAURA.jar"
GTF="/path/annotation.gtf"
OUTROOT="/path/APAURA_results"
LIBRARY="reverse"

run_comparison(){
  name="$1"; control="$2"; case="$3"; samples="$4"
  out="${OUTROOT}/${name}"
  mkdir -p "$out"

  java -jar "$APAURA" short-read --samples "$samples" --gtf "$GTF" --control "$control" --case "$case" --library "$LIBRARY" --out "$out"
}

run_comparison "comparison_1" "controlA" "caseA" "comparison_1.samples.tsv"
run_comparison "comparison_2" "controlB" "caseB" "comparison_2.samples.tsv"
```

Recommended practices:

- preserve every sample sheet with its result directory;
- use the same annotation across comparisons that will be compared;
- keep library-strand interpretation consistent for libraries generated with the same protocol;
- record the APAURA version;
- capture logs;
- avoid private absolute paths in scripts committed to GitHub.
