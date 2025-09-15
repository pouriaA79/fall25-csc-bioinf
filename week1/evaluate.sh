#!/usr/bin/env bash
set -euxo pipefail

# =========================
# =========================
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
WEEK1_DIR="$REPO_ROOT/week1"
CODE_DIR="$WEEK1_DIR/codes"      
DATA_DIR="$WEEK1_DIR/data"

PY_ENTRY="${PY_ENTRY:-$CODE_DIR/python/main.py}"
CODON_ENTRY="${CODON_ENTRY:-$CODE_DIR/codon/main_type.py}"

# Codon binary
if command -v codon >/dev/null 2>&1; then
  CODON_BIN="codon"
else
  CODON_BIN="${CODON_BIN:-$HOME/.codon/bin/codon}"
fi

CODON_ARGS_DEFAULT=(${CODON_ARGS_DEFAULT:-batch 25 100 0 1 3000 0 0 200000})

RESULTS_CSV="${RESULTS_CSV:-$WEEK1_DIR/results.csv}"
: > "$RESULTS_CSV"

# =========================
# =========================
fmt_duration () { local s="$1"; printf "%d:%02d:%02d" $((s/3600)) $(((s%3600)/60)) $((s%60)); }

run_and_time () {
  local __out_var="$1"; shift
  local start end; start="$(date +%s)"; "$@" || true; end="$(date +%s)"
  printf -v "$__out_var" "%s" "$(( end - start ))"
}

n50_of_fasta () {
  local fasta="$1"
  if [[ ! -s "$fasta" ]]; then echo "NA"; return 0; fi
  local lengths total n50
  lengths="$(awk '
    BEGIN{l=0}
    /^>/{ if(l>0) print l; l=0; next }
    { l += length($0) }
    END{ if(l>0) print l }
  ' "$fasta" || true)"
  [[ -z "$lengths" ]] && { echo "NA"; return 0; }
  total="$(echo "$lengths" | awk '{s+=$1} END{print s+0}')"
  [[ "$total" -eq 0 ]] && { echo "NA"; return 0; }
  n50="$(echo "$lengths" | sort -nr | awk -v T="$total" 'BEGIN{c=0; t=T/2}{c+=$1; if(c>=t){print $1; exit}}')"
  [[ -z "$n50" ]] && echo "NA" || echo "$n50"
}

find_contigs_fasta () {
  local dataset="$1"; local dspath="$2"; local out=""
  if   [[ -s "$dspath/contigs.fasta"        ]]; then out="$dspath/contigs.fasta"
  elif [[ -s "$dspath/contig_python.fasta"  ]]; then out="$dspath/contig_python.fasta"
  elif compgen -G "$dspath/contig_*.fasta" > /dev/null; then
    local merged="$dspath/merged_contigs.fasta"; rm -f "$merged"
    cat "$dspath"/contig_*.fasta > "$merged" || true
    [[ -s "$merged" ]] && out="$merged"
  fi
  if [[ -z "$out" ]]; then
    if   [[ -s "./${dataset}_batched/contigs.fasta" ]]; then out="./${dataset}_batched/contigs.fasta"
    elif [[ -s "./${dataset}_single/contigs.fasta"  ]]; then out="./${dataset}_single/contigs.fasta"
    elif [[ -s "./${dataset}/contigs.fasta"         ]]; then out="./${dataset}/contigs.fasta"
    else
      local merged="./${dataset}_merged_contigs.fasta"; rm -f "$merged"
      if compgen -G "./${dataset}/contig_*.fasta" > /dev/null; then
        cat "./${dataset}/contig_"*.fasta > "$merged"
      elif compgen -G "./${dataset}_*/contig_*.fasta" > /dev/null; then
        cat "./${dataset}_"/contig_*.fasta > "$merged" || true
      fi
      [[ -s "$merged" ]] && out="$merged"
    fi
  fi
  echo "$out"
}

discover_datasets () {
  local ds; for ds in data1 data2 data3; do
    [[ -d "$DATA_DIR/$ds" ]] && echo "$DATA_DIR/$ds"
  done
}

# =========================
# =========================
cd "$WEEK1_DIR"

printf "%-8s\t%-8s\t%-8s\t%-6s\n" "Dataset" "Language" "Runtime" "N50"
printf -- "-------------------------------------------------------------------\n"
echo "dataset,language,runtime_sec,runtime_hms,n50" >> "$RESULTS_CSV"

DATASETS=($(discover_datasets))
[[ ${#DATASETS[@]} -eq 0 && -d "$DATA_DIR/data1" ]] && DATASETS+=("$DATA_DIR/data1")

for ds_path in "${DATASETS[@]}"; do
  ds_name="$(basename "$ds_path")"

  # Python
  py_runtime=0
  if [[ -f "$PY_ENTRY" ]]; then
    run_and_time py_runtime python3 -u "$PY_ENTRY" "$ds_path"
    py_fa="$(find_contigs_fasta "$ds_name" "$ds_path")"
    py_n50="$(n50_of_fasta "$py_fa")"
    printf "%-8s\t%-8s\t%-8s\t%-6s\n" "$ds_name" "python" "$(fmt_duration "$py_runtime")" "$py_n50"
    echo "$ds_name,python,$py_runtime,$(fmt_duration "$py_runtime"),$py_n50" >> "$RESULTS_CSV"
  else
    printf "%-8s\t%-8s\t%-8s\t%-6s\n" "$ds_name" "python" "0:00:00" "NA"
    echo "$ds_name,python,0,0:00:00,NA" >> "$RESULTS_CSV"
  fi

  # Codon
  codon_runtime=0
  if [[ -f "$CODON_ENTRY" ]]; then
    run_and_time codon_runtime \
      "$CODON_BIN" run -release -plugin seq "$CODON_ENTRY" "$ds_path"
    codon_fa="$(find_contigs_fasta "$ds_name" "$ds_path")"
    codon_n50="$(n50_of_fasta "$codon_fa")"
    printf "%-8s\t%-8s\t%-8s\t%-6s\n" "$ds_name" "codon" "$(fmt_duration "$codon_runtime")" "$codon_n50"
    echo "$ds_name,codon,$codon_runtime,$(fmt_duration "$codon_runtime"),$codon_n50" >> "$RESULTS_CSV"
  else
    printf "%-8s\t%-8s\t%-8s\t%-6s\n" "$ds_name" "codon" "0:00:00" "NA"
    echo "$ds_name,codon,0,0:00:00,NA" >> "$RESULTS_CSV"
  fi
done

echo "CSV saved to $RESULTS_CSV"
