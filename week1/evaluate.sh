#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WEEK1_DIR="$SCRIPT_DIR"
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
  # usage: run_and_time_logged <outvar_seconds> <logfile> -- <cmd> <args...>
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
  # طول همه‌ی خطوط توالی بین هدرها جمع می‌شود
  local lengths total
  lengths="$(awk 'BEGIN{l=0}/^>/{if(l>0)print l;l=0;next}{l+=length($0)}END{if(l>0)print l}' "$fa")"
  [[ -n "$lengths" ]] || { echo NA; return; }
  total="$(echo "$lengths" | awk '{s+=$1}END{print s+0}')"
  [[ "$total" -gt 0 ]] || { echo NA; return; }
  echo "$lengths" | sort -nr | awk -v T="$total" 'BEGIN{c=0; t=T/2}{c+=$1; if(c>=t){print $1; exit}}'
}

discover_datasets(){ for ds in data1 data2 data3; do [[ -d "$DATA_DIR/$ds" ]] && echo "$DATA_DIR/$ds"; done; }

find_contigs_in_ds(){  # فقط داخل فولدر دیتاست
  local d="$1"
  if   [[ -s "$d/contigs.fasta"       ]]; then echo "$d/contigs.fasta"
  elif [[ -s "$d/contig_python.fasta" ]]; then echo "$d/contig_python.fasta"
  elif compgen -G "$d/contig_*.fasta" >/dev/null; then
    local m="$d/merged_contigs.fasta"; rm -f "$m"; cat "$d"/contig_*.fasta > "$m" || true; [[ -s "$m" ]] && echo "$m"
  fi
}

echo -e "Dataset\tLanguage\tRuntime\tN50"
echo "------------------------------------------------------"

cd "$WEEK1_DIR"
datasets=($(discover_datasets))
[[ ${#datasets[@]} -gt 0 ]] || { echo "No datasets under $DATA_DIR"; exit 1; }

for ds in "${datasets[@]}"; do
  name="$(basename "$ds")"

  # ============= Python =============
  py_rt=0; py_log="$LOG_DIR/py_${name}.log"
  if [[ -f "$PY_ENTRY" ]] && run_and_time_logged py_rt "$py_log" -- python3 -u "$PY_ENTRY" "$ds"; then
    py_fa="$(grep -m1 '^WROTE:' "$py_log" | awk '{print $2}')"
    [[ -z "${py_fa:-}" ]] && py_fa="$(find_contigs_in_ds "$ds")"
    if [[ -n "${py_fa:-}" ]]; then
      echo "[PY/$name] file: $py_fa"
      ls -l "$py_fa" || true
      echo "[PY/$name] size(bytes): $(wc -c <"$py_fa" || echo 0)"
    else
      echo "[PY/$name] contigs not found"
    fi
    py_n50="$(n50_of_fasta "${py_fa:-/dev/null}")"
    echo -e "$name\tpython\t$(fmt_duration "$py_rt")\t$py_n50"
    echo "$name,python,$py_rt,$(fmt_duration "$py_rt"),$py_n50" >> "$RESULTS_CSV"
  else
    echo "[PY/$name] FAILED (see $py_log)"
    echo "$name,python,0,0:00:00,NA" >> "$RESULTS_CSV"
  fi

  # ============= Codon =============
  codon_rt=0; codon_log="$LOG_DIR/codon_${name}.log"
  if [[ -f "$CODON_ENTRY" ]] && run_and_time_logged codon_rt "$codon_log" -- "$CODON_BIN" run -release -plugin seq "$CODON_ENTRY" "$ds"; then
    co_fa="$(grep -m1 '^WROTE:' "$codon_log" | awk '{print $2}')"
    [[ -z "${co_fa:-}" ]] && co_fa="$(find_contigs_in_ds "$ds")"
    if [[ -n "${co_fa:-}" ]]; then
      echo "[CODON/$name] file: $co_fa"
      ls -l "$co_fa" || true
      echo "[CODON/$name] size(bytes): $(wc -c <"$co_fa" || echo 0)"
    else
      echo "[CODON/$name] contigs not found"
    fi
    co_n50="$(n50_of_fasta "${co_fa:-/dev/null}")"
    echo -e "$name\tcodon\t$(fmt_duration "$codon_rt")\t$co_n50"
    # ✅ این‌جا قبلاً اشتباه "python" نوشته شده بود
    echo "$name,codon,$codon_rt,$(fmt_duration "$codon_rt"),$co_n50" >> "$RESULTS_CSV"
  else
    echo "[CODON/$name] FAILED (see $codon_log)"
    echo "$name,codon,0,0:00:00,NA" >> "$RESULTS_CSV"
  fi
done

echo "CSV saved to $RESULTS_CSV"
echo "Logs at: $LOG_DIR"
