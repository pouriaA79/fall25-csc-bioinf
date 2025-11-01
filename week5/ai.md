### AI TOOL AND MODEL

**ChatGPT**

---

### WHAT I USED AI FOR

During this project, I used AI mainly to **clarify concepts**, **debug pipeline issues**, and **understand tool behaviors** while keeping full control of the actual implementation.  
The main cases where I used it were:

1. **Clarifying genomic analysis steps**  
   I asked conceptual questions about how each stage of the pipeline works — for example, the logic behind *minimap2 alignment*, *bcftools variant calling*, and *HapCUT2 phasing*.  
   These explanations helped me understand what each parameter does and how read-based phasing connects to downstream star-allele interpretation.

2. **Understanding and automating IGV snapshots**  
   I spent considerable time getting **IGV headless snapshots** to work inside the notebook and CI environment.  
   The AI helped me understand how to write a proper IGV `.batch` script, set paths correctly, and use `Xvfb` so IGV could run without an interactive GUI.  
   This part was one of the trickiest parts of the assignment, and AI guidance was essential for solving the automation issues.

3. **Debugging CI errors and environment setup**  
   I asked for help when GitHub Actions failed (for example, with `libpython.so` errors, missing dependencies, or Codon memory limits).  
   The AI helped me fix those issues safely and optimize the workflow so the notebook executes end-to-end automatically.

---

### PROMPTS

1. *Why does HapCUT2 produce a 12-column block file, and how can I fix it so Whatshap can convert it to a phased VCF?*  
2. *How do I configure IGV to take automatic snapshots of discordant loci without a GUI?*  
3. *How should I interpret the phasing results (0|1 vs 1|0) when comparing to PharmVar star-alleles?*

---

### SUMMARY

I used AI as a **learning and troubleshooting assistant**. 
It helped me understand the bioinformatics logic behind each step, debug complex environment issues, and successfully automate IGV snapshots and CI runs.
