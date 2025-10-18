#!/usr/bin/env bash
set -e

echo "Method              Language    Runtime"
echo "--------------------------------------"

# مسیر اصلی پروژه
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

# مسیر داده‌ها
DATA_PATH="${DIR}/data"

if [ ! -d "$DATA_PATH" ]; then
    echo "❌ Error: data folder not found at $DATA_PATH"
    exit 1
fi

# --------------------------
# اجرای تست‌های Python
# --------------------------
echo "▶ Running Python tests..."
PY_START=$(date +%s%3N)

# اجرای مستقیم با نمایش خروجی در ترمینال
python3 tests/run_tests.py "$DATA_PATH"

PY_END=$(date +%s%3N)
PY_TOTAL=$((PY_END - PY_START))

echo ""
echo "python total runtime: ${PY_TOTAL}ms"
echo "--------------------------------------"

# --------------------------
# اجرای تست‌های Codon
# --------------------------
echo "▶ Running Codon tests..."
CODON_START=$(date +%s%3N)

# اجرای مستقیم با نمایش خروجی در ترمینال
codon run tests/run_tests_codon.py "$DATA_PATH"

CODON_END=$(date +%s%3N)
CODON_TOTAL=$((CODON_END - CODON_START))

echo ""
echo "codon total runtime: ${CODON_TOTAL}ms"
echo "--------------------------------------"

# --------------------------
# خلاصه نهایی
# --------------------------
echo ""
echo "Language    Total Runtime"
echo "--------------------------"
printf "python      %dms\n" "$PY_TOTAL"
printf "codon       %dms\n" "$CODON_TOTAL"
