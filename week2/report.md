# Week 2 Report — `trviz` Port (Python → Codon)

## Overview

The goal of this week's assignment was to select an alternative to the main BioPython task and port its logic to Codon. I chose the `trviz` alternative, which involves implementing the core non-visual components of the tandem repeat visualization pipeline.

**I successfully:**
* Implemented the `trviz` data processing pipeline (Decomposer, Encoder, Aligner) in a single, Codon-compatible script located in the `code/` folder.
* Fully implemented the pipeline which reads input data from a FASTA file and produces the final output as an image using the `Visualizer` component (a sample output image has been uploaded to the repository).
* Utilized Codon's Python bridge to interface with `matplotlib` for the visualization component, as permitted by the assignment.
* Developed a robust, self-contained test suite (`test.py`) that verifies the functionality of all core components within the Codon environment.
* Configured a GitHub Actions CI workflow (`week2-ci.yml`) that automatically installs all dependencies (including `mafft`) and runs the test suite on every push.

The project is fully functional and passes all tests in the CI environment.

---

## Porting Notes & Technical Challenges

The process of converting the `trviz` logic to a stable, testable Codon program involved several specific challenges related to compiler behavior, external tool interaction, and module resolution.

### 1. Codon/Matplotlib: `Image size too large` Error
* **Problem:** When processing FASTA files with long tandem repeat sequences, the `matplotlib` backend raised a `PyError` because the calculated image dimensions in pixels exceeded its internal limit (2^16).
* **Fix:** The solution involved two adjustments in the `Visualizer.trplot` function: reducing the width multiplier in the `figsize` calculation and lowering the save resolution from `dpi=300` to `dpi=150`.

### 2. Managing Static Typing in Codon
* **Problem:** As a compiled language, Codon requires explicit type hints for variables and function outputs. A key challenge during the porting process was identifying the precise types for some of the complex outputs, especially when interacting with Python libraries.
* **Fix:** This was resolved through careful inspection of outputs and the use of appropriate type hints from the `typing` module (e.g., `Tuple`, `Dict`, `List`) to make the code fully understandable to the Codon compiler.

### 3. Codon Compiler: `ModuleNotFoundError`
* **Problem:** This was the most persistent challenge. With the project structured into separate `code/` and `test/` directories, the Codon compiler consistently failed to find the `main_type` module from `test.py`.
* **Root Cause & Final Solution:** These failures indicated that Codon's compiler resolves module paths *before* runtime, ignoring any dynamic changes to `sys.path`. To completely bypass this compiler-specific issue, all necessary classes and functions were **copied directly into the `test.py` script**. This transformed the test script into a fully standalone program, definitively solving the import error.

---


## Automation (CI & Testing)

* **Testing (`test.py`):** A self-contained, Codon-only test suite was developed. It includes unit tests for `Decomposer` and the encoding functions, as well as a robust integration test for `MotifAligner` that verifies its interaction with `mafft`.
  **A key aspect of the testing strategy was the use of synthetic data.** Since no specific test dataset was provided for this task, simple, predictable data was created directly within the `test.py` script. This allowed for the isolated and precise testing of each component: unit tests for the `Decomposer` and encoding functions, and a robust integration test for the `MotifAligner`. All formal tests executed by the CI are performed on this synthetic data. 

* **CI (`week2-ci.yml`):** A GitHub Actions workflow automates the entire testing process. On every `push` and `pull_request` that affects the `week2/` directory, the CI job:
    1.  Checks out the repository code.
    2.  Installs the exact required version of Codon and the `seq` plugin.
    3.  Verifies the Codon installation with a version check and a smoke test.
    4.  Installs system dependencies (`mafft`) and Python libraries (`matplotlib`, `biopython`).
    5.  Configures the Codon-Python bridge.
    6.  Executes the standalone `week2/test/test.py` script.

This CI provides a complete, automated guarantee of the project's correctness.

---

## Time Estimate

* **Estimated Time to Complete:** (Approximately 30 hours)

