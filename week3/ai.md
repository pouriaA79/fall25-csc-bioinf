### AI TOOL AND MODEL

**Gemini — 2.5 Flash**

---

### WHAT I USED AI FOR

1. **Cleaning and optimizing the CI pipeline** with automatic Codon installation, Python bridge setup, and runtime benchmarking.
2. **Explaining differences between Python and Codon module resolution**, and fixing `ModuleNotFoundError` by correctly setting `CODONPATH`.
3. **Debugging timing logic** to ensure `evaluate.sh` outputs exact total runtimes for both Python and Codon tests.

---

### PROMPTS

1. *"How should I set up `CODONPATH` correctly?"*
2. *"How can I safely measure execution time in milliseconds inside Codon using the `time` module?"*
3. *"Explain how to mirror Python’s `sys.path` behavior in Codon when importing from sibling directories."*
4. *"Clean and simplify this GitHub Actions workflow that installs Codon, runs Python and Codon tests, and uploads benchmark results as artifacts."*

