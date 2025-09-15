#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
WEEK1_DIR="$REPO_ROOT/week1"
CODE_DIR="$WEEK1_DIR/codes"
DATA_DIR="$WEEK1_DIR/data"
LOG_DIR="$WEEK1_DIR/ci_logs"
mkdir -p "$LOG_DIR"

PY_ENTRY="${PY_ENTRY:-$CODE_DIR/python/main.py}"
CODON_ENTRY="${CODON_ENTRY:-$CODE_DIR/codon/main_type.py}"
CODON_BIN="$(command -v codon || echo "$HOME/.codon/bin/codon")"

RESULTS_CSV="$WEEK1_DIR/results.csv"
echo "dataset,language,runtime_sec,runtime_hms,n50" > "$RESULTS_CSV"

fmt_duration(){ s="$1"; printf "%d:%02d:%02d" $((s/3600)) $(((s%3600)/60)) $((s%60)); }

run_and_time_logged () {
  local __out="$1"; shift
  local log="$1"; shift
  [[ "${1:-}" == "--" ]] && shift || true
  local start end rc
  start="$(date +%s)"; set +e; "$@" >"$log" 2>&1; rc=$?; set -e; end="$(date +%s)"
  printf -v "$__out" "%s" "$((end-start))"
  return $rc
}

n50_of_fasta () {
  local fa="$1"; [[ -s "$fa" ]] || { echo NA; return; }
  local lengths total
  lengths="$(awk 'BEGIN{l=0}/^>/{if(l>0)print l;l=0;next}{l+=length($0)}END{if(l>0)print l}' "$fa")"
  [[ -n "$lengths" ]] || { echo NA; return; }
  total="$(echo "$lengths" | awk '{s+=$1}END{print s+0}')"
  [[ "$total" -gt 0 ]] || { echo NA; return; }
  echo "$lengths" | sort -nr | awk -v T="$total" 'BEGIN{c=0; t=T/2}{c+=$1; if(c>=t){print $1; exit}}'
}

discover_datasets(){ for ds in data1 data2 data3; do [[ -d "$DATA_DIR/$ds" ]] && echo "$DATA_DIR/$ds"; done; }

find_contigs_in_ds(){ # فقط داخل فولدر دیتاست
  local d="$1"
  if   [[ -s "$d/contigs.fasta"       ]]; then echo "$d/contigs.fasta"
  elif [[ -s "$d/contig_python.fasta" ]]; then echo "$d/contig_python.fasta"
  elif compgen -G "$d/contig_*.fasta" >/dev/null; then
    local m="$d/merged_contigs.fasta"; rm -f "$m"; cat "$d"/contig_*.fasta > "$m" || true; [[ -s "$m" ]] && echo "$m"
  fi
}

cd "$WEEK1_DIR"
datasets=($(discover_datasets))
[[ ${#datasets[@]} -gt 0 ]] || { echo "No datasets under $DATA_DIR"; exit 1; }

printf "%-8s\t%-8s\t%-8s\t%-6s\n" "Dataset" "Language" "Runtime" "N50"
printf -- "-------------------------------------------------------------------\n"

for ds in "${datasets[@]}"; do
  name="$(basename "$ds")"
  # پاک‌سازی خروجی‌های قدیمی
  rm -f "$ds/contigs.fasta" "$ds/merged_contigs.fasta" "$ds"/contig_*.fasta || true

  # -------- Python --------
  py_rt=0; py_log="$LOG_DIR/py_${name}.log"
  if [[ -f "$PY_ENTRY" ]] && run_and_time_logged py_rt "$py_log" -- python3 -u "$PY_ENTRY" "$ds"; then
    py_fa_log="$(grep -m1 '^WROTE:' "$py_log" | awk '{print $2}')"
    py_fa="${py_fa_log:-$(find_contigs_in_ds "$ds")}"
    echo "PY found: ${py_fa:-<none>}"
    py_n50="$(n50_of_fasta "$py_fa")"
    printf "%-8s\t%-8s\t%-8s\t%-6s\n" "$name" "python" "$(fmt_duration "$py_rt")" "$py_n50"
    echo "$name,python,$py_rt,$(fmt_duration "$py_rt"),$py_n50" >> "$RESULTS_CSV"
  else
    echo "Python FAILED on $name (see $py_log)"
    echo "$name,python,0,0:00:00,NA" >> "$RESULTS_CSV"
  fi

  # -------- Codon --------
  codon_rt=0; codon_log="$LOG_DIR/codon_${name}.log"
  if [[ -f "$CODON_ENTRY" ]] && run_and_time_logged codon_rt "$codon_log" -- "$CODON_BIN" run -release -plugin seq "$CODON_ENTRY" "$ds"; then
    co_fa_log="$(grep -m1 '^WROTE:' "$codon_log" | awk '{print $2}')"
    co_fa="${co_fa_log:-$(find_contigs_in_ds "$ds")}"
    echo "CODON found: ${co_fa:-<none>}"
    co_n50="$(n50_of_fasta "$co_fa")"
    printf "%-8s\t%-8s\t%-8s\t%-6s\n" "$name" "codon" "$(fmt_duration "$codon_rt")" "$co_n50"
    echo "$name,codon,$codon_rt,$(fmt_duration "$codon_rt"),$co_n50" >> "$RESULTS_CSV"
  else
    echo "Codon FAILED on $name (see $codon_log)"
    echo "$name,codon,0,0:00:00,NA" >> "$RESULTS_CSV"
  fi
done

echo "CSV saved to $RESULTS_CSV"
echo "Logs at: $LOG_DIR"
