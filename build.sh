#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

echo "Building APAURA..."

rm -rf build
mkdir -p build

# ============================================================
# Standard build: Maven
# ============================================================

if command -v mvn >/dev/null 2>&1; then
    mvn clean package

    if [[ ! -f target/APAURA.jar ]]; then
        echo "ERROR: Maven completed but target/APAURA.jar was not found." >&2
        exit 2
    fi

    cp target/APAURA.jar build/APAURA.jar

# ============================================================
# Fallback build: javac + jar
# ============================================================

else
    command -v javac >/dev/null 2>&1 || {
        echo "ERROR: Java/JDK is required." >&2
        exit 2
    }

    command -v jar >/dev/null 2>&1 || {
        echo "ERROR: jar command is required." >&2
        exit 2
    }

    mkdir -p build/classes

    find src/main/java -name '*.java' | sort > build/sources.txt

    if javac --help 2>&1 | grep -q -- '--release'; then
        javac --release 8 -d build/classes @build/sources.txt
    else
        javac -source 8 -target 8  -d build/classes  @build/sources.txt
    fi

    # Preserve packaged resources.
    if [[ -d src/main/resources ]]; then
        cp -R src/main/resources/. build/classes/
    fi

    jar cfe build/APAURA.jar org.apaura.APAURA -C build/classes .
fi

# ============================================================
# Create APAURA command-line launcher
# ============================================================

cat > apaura <<'LAUNCHER'
#!/usr/bin/env bash
set -euo pipefail

APAURA_HOME="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -f "$APAURA_HOME/build/APAURA.jar" ]]; then
    APAURA_JAR="$APAURA_HOME/build/APAURA.jar"
elif [[ -f "$APAURA_HOME/target/APAURA.jar" ]]; then
    APAURA_JAR="$APAURA_HOME/target/APAURA.jar"
else
    echo "ERROR: APAURA.jar was not found." >&2
    echo "Run ./build.sh first." >&2
    exit 2
fi

if [[ -n "${JAVA_HOME:-}" && -x "$JAVA_HOME/bin/java" ]]; then
    JAVA="$JAVA_HOME/bin/java"
elif command -v java >/dev/null 2>&1; then
    JAVA="$(command -v java)"
else
    echo "ERROR: Java runtime was not found." >&2
    exit 2
fi

exec "$JAVA" -jar "$APAURA_JAR" "$@"
LAUNCHER

chmod +x apaura

# Remove the old root-level JAR, if present.
rm -f APAURA.jar

echo
echo "============================================"
echo " APAURA build complete"
echo "============================================"
echo
echo "Executable:"
echo "  $ROOT/apaura"
echo
echo "Examples:"
echo "  ./apaura --help"
echo "  ./apaura run-3seq --help"
echo "  ./apaura short-read --help"
echo
