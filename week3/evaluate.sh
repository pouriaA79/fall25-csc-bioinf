
set -e

echo "Language    Total Runtime"
echo "--------------------------"

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

extract_total_runtime() {
    total=0
    while IFS= read -r line; do
        if [[ "$line" =~ Runtime:\ ([0-9]+)ms ]]; then
            total=$((total + BASH_REMATCH[1]))
        fi
    done
    echo $total
}
export PYTHONPATH="$DIR/code"
PY_OUTPUT=$(python3 tests/test_distances_python.py; \
            python3 tests/test_neighbor_joining_python.py; \
            python3 tests/test_upgma_python.py)

PY_TOTAL=$(echo "$PY_OUTPUT" | extract_total_runtime)
echo "python      ${PY_TOTAL}ms"

unset PYTHONPATH

export CODONPATH="$DIR/code"
CODON_OUTPUT=$(codon run tests/test_distances.codon; \
               codon run tests/test_neighbor_joining.codon; \
               codon run tests/test_upgma.codon)

CODON_TOTAL=$(echo "$CODON_OUTPUT" | extract_total_runtime)
echo "codon       ${CODON_TOTAL}ms"


