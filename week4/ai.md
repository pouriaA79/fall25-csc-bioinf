
### AI TOOL AND MODEL

**Femini — Flash-2.5**

---

### WHAT I USED AI FOR

1. **Clarifying alignment algorithm logic**, especially for *Affine-gap* and *Semi-global alignment*, where I asked conceptual questions to better understand how scoring matrices and gap penalties are computed.
2. **Debugging CI errors** in GitHub Actions — I shared the output of the failed runs (like the `libpython.so` and memory-limit errors) and asked how to fix them safely.
3. **Improving the workflow configuration**, particularly ensuring Codon and Python are correctly linked and dependencies are automatically detected inside the CI runner.

---

### PROMPTS

1. *This CI run fails with `libpython.so: cannot open shared object file`. What should I do to fix it?*
2. *Can you explain in detail how affine-gap alignment matrices (M, Ix, Iy) are updated each step?*
3. *Why does my GitHub Actions runner abort Codon runs without an error, and how can I reduce memory usage?*

