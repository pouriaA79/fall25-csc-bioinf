# Week 1 Report — de Bruijn Assembler (Python → Codon)

## Overview

The goal of this week was to (1) set up automation with GitHub Actions, (2) port a simple de Bruijn graph genome assembler from Python to Codon and run it, and (3) reproduce results (N50) on the provided datasets.

I successfully:

* Ran the Python baseline on **data1–data3**
* Ported the code to **Codon** and matched **N50** values with Python on **data1–data3**
* Automated both runs with **CI** and produced a `results.csv` plus contig artifacts

`data4` currently **fails** for both Python and Codon under my environment (explained below).

---
## Porting Notes (Python → Codon) & Fixes

1. **Unnecessary matplotlib import**

   * Early runs failed due to a stray `matplotlib` import that wasn’t used.
     **Fix:** removed/commented it.

2. **OS functions in Codon**

   * Some `os` features aren’t implemented natively in Codon.
   * **Fix:** use the Python bridge explicitly:

     ```python
     from python import os, sys
     ```

     so `os.path` and friends behave like CPython.

3. **Type realization errors**

   * Codon needs concrete types during realization. In `dbg.py`, class fields had no explicit types, which caused errors like:

     ```
     dict.codon: error: 'NoneType' object has no attribute '__hash__'
     ```
   * **Fix:** added explicit type annotations, e.g.:

     ```python
     class Node:
         _children: Set[int]
         _count: int
         kmer: str
         visited: bool
         depth: int
         max_depth_child: Optional[int]
     ```
   * I initially tried mixing Python decorators (`@python`) to wrap ambiguous code paths, but the cleanest, most stable solution was to **fully annotate** the data structures.

## Automation (`evaluate.sh`)

* Runs both Python and Codon for `data1..data3`
* Finds the best output FASTA inside each dataset folder
* Computes **N50** (implemented in a small Python one-liner inside the script for robustness)
* Produces a CSV with columns:

  ```
  dataset,language,runtime_sec,runtime_hms,n50
  ```
* Uploads results and contigs as artifacts

---

## Results

* For **data1–data3**, **N50 matches** between Python and Codon (as intended).
* **Runtime** with Codon is **faster** (approximately \~½ of Python on my runs; see `week1/results.csv` in CI artifacts for the exact numbers captured on the runner).

> The exact values are available in the CI artifact `week1-results` (`week1/results.csv`).
> The `contigs` artifact contains the FASTA outputs per dataset.

---

## Data4 (Failure Analysis)

* Both Python and Codon **crashed** on `data4` under my environment.
* Symptoms observed:

  * Recursion depth / stack exhaustion on Python
  * WSL resource exhaustion / crash when increasing recursion limit or stack via `ulimit -s`
  * Large memory footprint for the de Bruijn graph on data4 scale
* Attempts:

  * Increased recursion limit (Python)
  * `ulimit -s` in WSL (stack size); WSL still crashed
* Conclusion:

  * With the current simple assembler and available resources, `data4` was **not reproducible** for me.
  * Per the assignment updates: I clearly document this in the report and mark the repo as not fully reproducible for `data4` as-is.

**Potential next steps (future work):**

* Iterative DFS everywhere (no recursion)
---

## Reproducibility Notes

* I reproduced N50 parity on **data1–data3** (Python vs. Codon).
* I clearly document the failure on **data4**.

