# **AI TOOL AND MODEL**

**ChatGPT**

---

# **WHAT I USED AI FOR**

I used AI only for **clarifying concepts** and **debugging technical issues** in the single-cell pipeline.
Specifically:

1. **Understanding pipeline steps**
   Getting explanations for how simpleaf / alevin-fry work and how their outputs map to Scanpy/AnnData.

2. **Fixing CI path and runtime errors**
   Debugging missing FASTQ files, incorrect directories, whitelist issues, and locating `quants_mat.mtx` in GitHub Actions.

3. **Improving reproducibility**
   Ensuring automatic execution, stable file paths, and correct artifact outputs.

All analysis, coding, and decisions were done by me; AI only helped with troubleshooting.

---

# **PROMPTS**

Examples of prompts I used:

* “Why is alevin-fry not finding `quants_mat.mtx` in CI?”
* “How should I pass whitelist and R1/R2 paths to simpleaf quant?”
* “How do I fix gdown permission errors in GitHub Actions?”

---

# **SUMMARY**

AI served as a **debugging and explanation assistant**.
I implemented the full workflow myself; AI only helped me understand tool behavior and resolve technical issues.
