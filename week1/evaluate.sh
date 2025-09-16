#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WEEK1_DIR="$SCRIPT_DIR"
CODE_DIR="$WEEK1_DIR/codes"
DATA_DIR="$WEEK1_DIR/data"

PY_ENTRY="${PY_ENTRY:-$CODE_DIR/python/main.py}"
CODON_ENTRY="${CODON_ENTRY:-$CODE_DIR/codon/main_type.py}"
CODON_BIN="$(command -v codon || echo "$HOME/.codon/bin/codon")"

RESULTS_CSV="$WEEK1_DIR/results.csv"
echo "dataset,language,runtime_sec,n50" > "$RESULTS_CSV"

fmt_duration(){ s="$1"; printf "%d:%02d:%02d" $((s/3600)) $(((s%3600)/60)) $((s%60)); }

run_and_time () {
  # usage: run_and_time <outvar_seconds> -- <cmd> <args...>
  local __out="$1"; shift
  [[ "${1:-}" == "--" ]] && shift || true
  local start end
  start="$(date +%s)"
  "$@"
  end="$(date +%s)"
  printf -v "$__out" "%s" "$((end-start))"
}

# پیدا کردن بهترین فایل FASTA داخل فولدر دیتاست
find_best_fasta_in_ds () {
  local d="$1"
  # 1) contigs.fasta
  [[ -s "$d/contigs.fasta" ]] && { echo "$d/contigs.fasta"; return 0; }
  # 2) contig_*.fasta → merge
  if compgen -G "$d/contig_*.fasta" >/dev/null; then
    local m="$d/merged_contigs.fasta"
    cat "$d"/contig_*.fasta > "$m" || true
    [[ -s "$m" ]] && { echo "$m"; return 0; }
  fi
  # 3) هر *.fasta موجود (بزرگ‌ترین)
  local best="" bestsz=0
  while IFS= read -r -d '' f; do
    local sz; sz=$(wc -c <"$f" || echo 0)
    if (( sz > bestsz )); then best="$f"; bestsz="$sz"; fi
  done < <(find "$d" -maxdepth 1 -type f -name "*.fasta" -print0 | sort -z)
  [[ -n "$best" ]] && { echo "$best"; return 0; }
  echo ""
}

# محاسبه N50 با پایتون (robust)
n50_of_fasta () {
  local fa="$1"
  [[ -s "$fa" ]] || { echo NA; return; }
  python3 - "$fa" <<'PY'
import sys
from pathlib import Path
fa = Path(sys.argv[1])
L = []
if not fa.is_file():
    print("NA"); sys.exit(0)
l = 0
with fa.open() as f:
    for line in f:
        if line.startswith('>'):
            if l>0:
                L.append(l); l=0
        else:
            l += len(line.strip())
if l>0: L.append(l)
if not L:
    print("NA"); sys.exit(0)
total = sum(L); half = total/2
s = 0
for x in sorted(L, reverse=True):
    s += x
    if s >= half:
        print(x); sys.exit(0)
print("NA")
PY
}

# دیتاست‌ها: data1..data3
datasets=()
for ds in data1 data2 data3; do
  [[ -d "$DATA_DIR/$ds" ]] && datasets+=("$DATA_DIR/$ds")
done
[[ ${#datasets[@]} -gt 0 ]] || { echo "No datasets under $DATA_DIR"; exit 1; }

echo -e "Dataset\tLanguage\tRuntime\tN50"
echo "----------------------------------------------"

for ds in "${datasets[@]}"; do
  name="$(basename "$ds")"

  # ---------- Python ----------
  py_rt=0
  if [[ -f "$PY_ENTRY" ]]; then
    run_and_time py_rt -- python3 -u "$PY_ENTRY" "$ds"
    py_fa="$(find_best_fasta_in_ds "$ds")"
    py_n50="$(n50_of_fasta "${py_fa:-/dev/null}")"
    echo -e "$name\tpython\t$(fmt_duration "$py_rt")\t${py_n50:-NA}"
    echo "$ds_name,python,$py_runtime,$py_n50" >> "$RESULTS_CSV"
  else
    echo -e "$name\tpython\t0:00:00\tNA"
    echo "$name,python,0,0:00:00,NA" >> "$RESULTS_CSV"
  fi

  # ---------- Codon ----------
  co_rt=0
  if [[ -f "$CODON_ENTRY" ]]; then
    run_and_time co_rt -- "$CODON_BIN" run -release -plugin seq "$CODON_ENTRY" "$ds"
    co_fa="$(find_best_fasta_in_ds "$ds")"
    co_n50="$(n50_of_fasta "${co_fa:-/dev/null}")"
    echo -e "$name\tcodon\t$(fmt_duration "$co_rt")\t${co_n50:-NA}"
    echo "$ds_name,codon,$codon_runtime,$codon_n50" >> "$RESULTS_CSV"
  else
    echo -e "$name\tcodon\t0:00:00\tNA"
    echo "$name,codon,0,0:00:00,NA" >> "$RESULTS_CSV"
  fi
done

echo "CSV saved to $RESULTS_CSV"
