#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

# Prefer Maven when it is available.
if command -v mvn >/dev/null 2>&1; then
    mvn -q clean package
    cp target/APAURA.jar APAURA.jar
    echo "Built: APAURA.jar"
    exit 0
fi

# Otherwise build directly with the JDK.
command -v javac >/dev/null 2>&1 || {
    echo "ERROR: Java/JDK is required." >&2
    exit 2
}

command -v jar >/dev/null 2>&1 || {
    echo "ERROR: jar command is required." >&2
    exit 2
}

rm -rf build
mkdir -p build/classes

find src/main/java -name '*.java' | sort > build/sources.txt

if javac --help 2>&1 | grep -q -- '--release'; then
    javac --release 8 -d build/classes @build/sources.txt
else
    javac -source 8 -target 8 -d build/classes @build/sources.txt
fi

jar cfe APAURA.jar org.apaura.APAURA -C build/classes .

echo "Built: APAURA.jar"
