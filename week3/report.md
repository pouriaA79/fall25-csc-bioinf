# Week 3 Report — `biotite.phylo` Port (Python → Codon)

---

## **Overview**

The goal of this week’s assignment was to **analyze and port the `biotite.phylo` module** — a key part of the Biotite bioinformatics library — from **Python/Cython** to **Codon**, and then compare their runtime performances.

The `phylo` module implements algorithms for computing evolutionary **distances** and constructing **phylogenetic trees** using methods such as **UPGMA** and **Neighbor Joining**.

I successfully:

* **Ported all required algorithms** (`distance`, `upgma`, `neighbor_joining`) to Codon.
* **Created Codon-compatible test scripts** that replicate Biotite’s original Python tests.
* **Developed a unified benchmarking script** (`evaluate.sh`) that executes both Python and Codon tests and aggregates total runtimes in milliseconds.
* **Configured a GitHub Actions workflow** (`week3-ci.yml`) that installs dependencies, runs the tests automatically, and uploads performance artifacts.
* **Validated correctness and performance** — Codon produced identical results while achieving a significant speed improvement. Additionally, I verified the correctness of the algorithms using a small test example.

The final CI run passes all tests successfully and outputs the required runtime comparison table.

---

## **Porting Notes & Technical Challenges**

Porting from **Cython** to **Codon** required in-depth understanding of the original implementation, since Codon doesn’t directly support Cython’s optimized data structures.

### **1. Porting from Cython to Codon**

* **Problem:** Biotite’s original code used Cython-specific syntax, optimized loops, and typed memory views that are unsupported in Codon.
* **Fix:** Rewrote all functions using standard Codon types (`list`, `tuple`, and `dict`) and reimplemented tree logic using custom classes for nodes and edges.
  This made the code portable and readable while preserving computational correctness.

---

### **2. Import & Path Resolution Errors**

* **Problem:** Codon failed to locate the local `phylo` module during compilation, showing:

  ```
  error: no module named 'phylo'
  ```
* **Fix:** I used
  ```
  from .. code.phylo import upgma, Tree
  ```

### **3. Timing Measurement Consistency**

* **Problem:** Runtime differences were too small to observe on minimal input data.
* **Fix:** Increased dataset size and repetition count so that measured execution times became stable and statistically meaningful.
  Each test internally measures and reports its execution time using:

  ```python
  start = time.perf_counter()
  ...
  print(f"Runtime: {int((time.perf_counter() - start) * 1000)}ms")
  ```

---

### **4. CI Environment Setup**

* **Problem:** The GitHub Actions runner initially lacked Codon and Biotite installations.
* **Fix:** Updated the CI workflow to automatically:

  * Install **Codon v0.19.3** and the **seq plugin**
  * Install **Python 3.11** and dependencies: `biotite`, `biopython`, and `find_libpython`

---

## **Automation (Testing & CI)**

### **Testing**

All test cases (`test_distances`, `test_upgma`, and `test_neighbor_joining`) were ported and adapted for Codon.
Each test prints its runtime, enabling direct comparison between Python and Codon.

### **Benchmarking (`evaluate.sh`)**

The evaluation script:

1. Runs all **Python** tests and extracts total runtime.
2. Runs all **Codon** tests and extracts total runtime.
3. Produces a final summary table:

   ```
   Language    Total Runtime
   --------------------------
   python      5982ms
   codon       1249ms
   ```

### **GitHub Actions Workflow (`week3-ci.yml`)**

The CI pipeline automatically:

1. Checks out the repository.
2. Installs Codon and the `seq` plugin.
3. Sets up Python and installs dependencies.
4. Runs `evaluate.sh` inside the `week3` directory.
5. Uploads results as a GitHub artifact (`week3-results`).

This ensures that every push or pull request is automatically verified for correctness and performance.

---

## **Performance Results**

**Final Benchmark Output:**

```
Language    Total Runtime
--------------------------
python      5982ms
codon       1249ms
```
This confirms Codon’s effectiveness for computationally intensive bioinformatics algorithms.

---

## **Time Estimate**

| Task                              | Duration      |
| --------------------------------- | ------------- |
| Analyzing Biotite source          | 8 hours       |
| Porting algorithms                | 17 hours       |
| Fixing imports and runtime issues | 8 hour        |
| CI setup and verification         | 5 hour        |
| **Total Estimated Time**          | **≈ 38 hours** |

---

## **Conclusion**

This project demonstrates that **Codon** can serve as a high-performance alternative to Python for scientific computing.
By porting the `biotite.phylo` module, we achieved a major **performance gain** while keeping the implementation clean, modular, and easy to maintain.

All requirements for **Week 3** were met successfully.
The Codon implementation is fast, correct, and fully integrated into an automated testing and benchmarking pipeline.
