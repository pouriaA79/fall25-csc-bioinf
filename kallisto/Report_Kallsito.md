## ABSTRACT

Kallisto is an RNA-seq quantification program. It pseudo-aligns reads to a reference, producing a list of transcripts that are compatible
with each read while avoiding alignment of individual bases. Given all the short reads from the RNA-seq, it tells us how many came from each transcript, and as a result estimates how much each transcript is expressed in a sample. Our project re-implements the core idea of Kallisto in the Codon language and compares its output to the original kallisto. Our implementation builds a transcript de Bruijn graph-style index, pseudo-aligns paired-end reads, forms equivalence classes, and runs an EM algorithm to estimate transcript counts and TPM (Transcripts Per Million). We then compare its output to the original kallisto on a small test dataset.
## INTRODUCTION

In experiments, we want to know which genes/transcripts are on and how much they are expressed. We start by sequencing RNA. After sequencing we get lots of short reads, each read is a tiny piece of some transcript. The goal is to determine from all these reads which transcripts (or genes) are more expressed in the sample. A straightforward solution is to use a full aligner that aligns every read to a reference genome or transcript with base-by-base alignment, this works but can be slow, and it also gives more information than we need for expression quantification. We don’t really need the exact alignment, we just need to know which transcripts are compatible with the read. Kallisto was introduced to solve this problem. It builds a transcript de Bruijn graph (tDBG) by breaking transcripts into k-mers and organizing them into a graph. Then For each read, looks at its k-mers, use the index to find which transcripts contain those k-mers, and output the set of compatible transcripts instead of an exact alignment, this is called pseudo-alignment. Then it groups reads with the same transcript set into equivalence classes. Finally, Kallisto runs an EM algorithm to split those counts among transcripts. We re-implemented the main kallisto pipeline in Codon. It follows the same idea: index, pseudo-alignment, equivalence classes, EM, TPM. Then ran both tools on the same dataset and compare results, looking for where they match and where they differ.

## MATERIALS AND METHODS

### Overview
There are four main steps in our kallisto-style pipeline. First, build the tDBG index from transcripts.fasta using dbg_index.py. Second load the index and paired-end reads into quant_dbg.py, pseudo-align each fragment by mapping its k-mers through the index and computing the set of compatible transcripts. Third is to group fragments that share the same transcript set into equivalence classes and count how many fragments fall into each class. Finally, feeding the equivalence class counts, transcript lengths, and fragment length into an EM algorithm to estimate transcript-level fragment counts and converts these counts and effective lengths into TPM values.

### Data
For this project we used the small toy dataset. It has a FASTA file with 14 transcripts, two FASTQ files with paired-end reads, each read pair is a noisy little window into one transcript. The goal is from these read pairs, we need to estimate how much each transcript is expressed. Both kallisto and our Codon code take the same input files and produce counts and TPM. In the end we compare the tables of est_counts and TPM that each method produces for these transcripts.

### Kallisto Baseline
Before running the Codon, we first ran the original kallisto on this dataset and used its output as a baseline.

First we built the index:
``` 
kallisto index -i out_kallisto/ref.idx data/transcripts.fasta 
```
This command reads the 14 transcript sequences from data/transcripts.fasta and builds kallisto’s internal transcript de Bruijn graph (tDBG) index, which it saves in out_kallisto/ref.idx. This index is basically a compact structure that lets kallisto quickly look up which transcripts contain a given k-mer when it pseudoaligns reads.

Then we ran quantification:
```kallisto quant \
  -i out_kallisto/ref.idx \
  -o out/kallisto_ref \
  data/reads_1.fastq data/reads_2.fastq
 ```
 Here kallisto uses the index and the two FASTQ files to pseudoalign each paired-end read, group fragments into equivalence classes, and run its EM algorithm to estimate expression. The main output we care about is out/kallisto_ref/abundance.tsv

### Codon Implementation
The Codon implementation includes three files as below:
- dbg_index.py: builds an index (a small transcript de Bruijn graph–like structure) from the transcript FASTA.
- quant_dbg.py: loads the index and reads, does pseudoalignment, builds equivalence classes, runs EM, and computes TPM.
- common.py: helper functions for reading FASTQ files, computing effective length, and writing the abundance.tsv file.

### 1) Transcript de Bruijn graph index (dbg_index.py)

To build the index in Codon, we run:
```
codon run dbg_index.py -- data/transcripts.fasta 31 out_dbg/kalli31

```

It reads all transcript sequences from data/transcripts.fasta,
breaks each transcript into overlapping k-mers of length 31,
records for each k-mer, which transcript IDs contain it,
groups k-mers into simple paths called “unitigs”, and
writes out three main mappings:
- k and the number of unitigs (.tdbg.meta)
- unitig_id → transcript IDs (.tdbg.unitigs.tsv)
- k-mer → unitig_id (.tdbg.k2u.tsv)
- plus a file with transcript lengths (.tids.txlen.tsv).

Later, quant_dbg.py uses these files to build a dictionary k2tix that maps each k-mer string directly to the indices of transcripts that contain it. This is the main structure used for pseudoalignment.

### 2) Pseudoalignment and equivalence classes (quant_dbg.py)

Then we run the Codon quantification step like this:
```
export CODON_ARGS="out_dbg/kalli31 100 \
  data/reads_1.fastq data/reads_2.fastq \
  out_dbg/dummy --mode union \
  --min-khits 1 --skip-budget 3 \
  --bootstrap 50 --seed 123"

codon run quant_dbg.py
 ```

Here we tell the script:
- which index to use (out_dbg/kalli31)
- the fragment length (100-150)
- the paired-end reads (reads_1.fastq, reads_2.fastq), and
that we want union mode, at least one k-mer hit, and 50 bootstrap replicates.

Inside quant_dbg.py, the script
loads the index and builds k2tix,
it reads all sequences from the FASTQ files using read_fastq from common.py. 
For each read, it slides a window of length k along the read, looks up each k-mer in k2tix, collects all transcript indices that share any of those k-mers. In union mode, merges them into one sorted list, this is the pseudoalignment for that read. For paired-end reads, it pseudoaligns each mate and takes the union of their compatible transcript sets,then groups fragments with the same sorted transcript set into an equivalence class and keeps a count of how many fragments fall into each class.

### 3) EM-based quantification and TPM
Given the equivalence classes, the script estimates transcript abundances with an EM algorithm (em_quant):
It first computes an effective length for each transcript, 
effective_length(L, frag_len) is basically max(1, L − frag_len + 1).
This corrects for the fact that longer transcripts naturally generate more fragments.
It initializes a vector of transcript abundances theta based on these effective lengths.
Then each EM iteration splits each class’s count across its member transcripts according to the current theta and effective lengths, sums these allocations for each transcript, then
normalizes to get an updated theta. This repeats until theta stabilizes.
After convergence, the script computes estimated fragment counts est_counts per transcript from the final theta and
converts them to TPM by dividing them by effective length. Finally, it writes one row per transcript to abundance.tsv with following columns: target_id, length, eff_length, est_counts, tpm. This has the same columns as kallisto’s abundance.tsv, so we can compare them directly.


## RESULT
### Overall comparison with kallisto
The ranking of transcripts by expression is very similar between our Codon version and kallisto. For example, ENST00000040584.5 and ENST00000282507.7 are among the top expressed transcripts in both Codon and kallisto. Several other top transcripts also appear in roughly the same order, which means our implementation is capturing the main signal in the data.
Even though the raw estimated counts (est_counts) can differ quite a bit, the TPM values are much closer. This is what matters most biologically: which transcripts are on and which ones express more.

For example:

- ENST00000513300.5
Codon: TPM ≈ 1.6k, 
kallisto: TPM ≈ 11k
- ENST00000504685.5
Codon: TPM ≈ 515, 
kallisto: TPM ≈ 10k


So our Codon version is giving these two a lot less expression than kallisto does, and that mass is basically going to other overlapping transcripts instead. That’s not surprising, these are isoforms that share many k-mers, so reads are highly ambiguous. Our pipeline uses a simpler model, while kallisto has a more detailed fragment length model and a more complex implementation. Because of that, when a read could belong to multiple very similar transcripts, we don’t split the counts in exactly the same way.

### Codon result

| target_id         | length | eff_length  | est_counts | tpm         |
| ----------------- | ------ | ----------- | ---------- | ---------|
| ENST00000243056.4 | 2423   | 2324.000000 | 19.000000  | 3433.59 |
| ENST00000243082.4 | 2039   | 1940.000000 | 28.000000  | 6061.60 |
| ENST00000303460.4 | 1936   | 1837.000000 | 26.000000  | 5944.23 |
| ENST00000312492.2 | 1805   | 1706.000000 | 115.000000 | 28310.68 |


### Kallisto result

| target_id         | length | eff_length | est_counts | tpm       |
| ----------------- | ------ | ---------- | ---------- | --------- |
| ENST00000243056.4 | 2423   | 2245.98    | 42         | 3553.05   |
| ENST00000243082.4 | 2039   | 1861.98    | 55         | 5612.36   |
| ENST00000303460.4 | 1936   | 1758.98    | 47         | 5076.85   |
| ENST00000312492.2 | 1805   | 1627.98    | 228        | 26609.9   |


## Conlusion 

In this project, we reimplemented a simplified kallisto-style RNA-seq quantification pipeline in Codon, from building a de Bruijn graph index to k-mer–based pseudoalignment, equivalence classes, and EM-based abundance estimation. Using a small dataset, we compared our Codon implementation directly against the official kallisto output.
Overall, the two methods agree well on the main expression signal. The ranking of highly expressed transcripts and the general TPM patterns are very similar, and many transcripts have comparable TPM values in both tools. The biggest difference appear in isoforms that share many k-mers, where reads are highly ambiguous. In these cases, kallisto and our pipeline distribute counts differently, which is expected given that our implementation uses a simplified model, no bias correction, and a less optimized pseudoalignment procedure.
Despite these simplifications, our reimplementation reproduces the core behavior of kallisto on this dataset and shows that the pseudoalignment framework can still yield reasonable expression estimates in a minimal, educational setting.
