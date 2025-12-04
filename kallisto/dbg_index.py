# dbg_index.py — Fixed k2u mapping (Codon-safe CLI)

import sys, os

def _cli_args() -> list[str]:
    try:
        from python import sys as pysys
        if hasattr(pysys, "argv") and pysys.argv:
            return [str(x) for x in pysys.argv]
    except:
        pass
    try:
        import sys as codonsys
        if hasattr(codonsys, "argv") and codonsys.argv:
            return [str(x) for x in codonsys.argv]
        if hasattr(codonsys, "args") and codonsys.args:
            return [str(x) for x in codonsys.args]
    except:
        pass
    try:
        raw = os.getenv("CODON_ARGS")
        if raw is not None and raw.strip() != "":
            return raw.strip().split()
    except:
        pass
    try:
        return sys.argv
    except:
        return []

def is_acgt(s: str) -> bool:
    for ch in s:
        if ch not in ("A","C","G","T"):
            return False
    return True

def read_fasta(path: str) -> dict[str, str]:
    seqs: dict[str, str] = {}
    tid = ""
    seq = ""
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if tid != "":
                    seqs[tid] = seq.upper().replace("U", "T")
                tid = line[1:].split()[0]
                seq = ""
            else:
                seq += line
    if tid != "":
        seqs[tid] = seq.upper().replace("U", "T")
    return seqs

def build_dbg(tx: dict[str,str], k: int):
    indeg: dict[str,int] = {}
    outdeg: dict[str,int] = {}
    k2tids: dict[str,list[str]] = {}
    for tid, seq in tx.items():
        n = len(seq)
        if n < k: continue
        for i in range(n-k+1):
            km = seq[i:i+k]
            if not is_acgt(km): continue
            pref, suff = km[:k-1], km[1:]
            outdeg[pref] = outdeg[pref] + 1 if pref in outdeg else 1
            indeg[suff] = indeg[suff] + 1 if suff in indeg else 1
            if km in k2tids:
                k2tids[km].append(tid)
            else:
                k2tids[km] = [tid]
    for kmer in list(indeg.keys()):
        if kmer not in outdeg: outdeg[kmer] = 0
    for kmer in list(outdeg.keys()):
        if kmer not in indeg: indeg[kmer] = 0
    return indeg, outdeg, k2tids

def unique_sorted(xs: list[str]) -> list[str]:
    if not xs: return []
    xs.sort()
    out = [xs[0]]
    for x in xs[1:]:
        if x != out[-1]:
            out.append(x)
    return out

def compact_to_unitigs(indeg, outdeg, k, kmer2tids):
    prefix2kmers: dict[str,list[str]] = {}
    for km in kmer2tids.keys():
        pref = km[:k-1]
        if pref in prefix2kmers:
            prefix2kmers[pref].append(km)
        else:
            prefix2kmers[pref] = [km]
    for p in prefix2kmers:
        prefix2kmers[p] = unique_sorted(prefix2kmers[p])

    visited: set[str] = set()
    unitigs: list[list[str]] = []  
    k2u: dict[str,int] = {}

    for km in kmer2tids.keys():
        if km in visited: 
            continue


        uid = len(unitigs)
        tids_acc: list[str] = []
        kmers_this_unitig: list[str] = []

        cur = km
        while True:
            visited.add(cur)
            kmers_this_unitig.append(cur)
            tids_acc.extend(kmer2tids[cur])
            next_candidates = prefix2kmers[cur[1:]] if cur[1:] in prefix2kmers else []
            if len(next_candidates) != 1:
                break
            nxt = next_candidates[0]
            if (indeg[nxt[:k-1]] if nxt[:k-1] in indeg else 0) != 1:
                break
            if nxt in visited:
                break
            cur = nxt

        tids_acc = unique_sorted(tids_acc)
        unitigs.append(tids_acc)
        for x in kmers_this_unitig:
            k2u[x] = uid

    return unitigs, k2u

def write_tx_lengths(prefix, tx):
    with open(prefix + ".tids.txlen.tsv", "w") as w:
        for tid, seq in tx.items():
            w.write(f"{tid}\t{len(seq)}\n")

def write_index(prefix, k, unitigs, k2u):
    with open(prefix + ".tdbg.meta", "w") as w:
        w.write(f"k={k}\tn_unitig={len(unitigs)}\n")
    with open(prefix + ".tdbg.unitigs.tsv", "w") as w:
        for i, tids in enumerate(unitigs):
            w.write(f"{i}\t{','.join(tids)}\n")
    with open(prefix + ".tdbg.k2u.tsv", "w") as w:
        for km, uid in k2u.items():
            w.write(f"{km}\t{uid}\n")

def main():
    raw = _cli_args()[1:]
    args = [a for a in raw if a and a != "--"]
    if len(args) < 3:
        print("Usage:\n  codon run dbg_index.py -- <transcripts.fasta> <k> <out_prefix>")
        print("Args seen:", args)
        sys.exit(1)

    fasta, k, outp = args[0], int(args[1]), args[2]
    tx = read_fasta(fasta)
    indeg, outdeg, k2tids = build_dbg(tx, k)
    unitigs, k2u = compact_to_unitigs(indeg, outdeg, k, k2tids)
    write_index(outp, k, unitigs, k2u)
    write_tx_lengths(outp, tx)
    print(f"[dbg-index] k={k}, unitigs={len(unitigs)}, kmers={len(k2u)}")
    print(f"[dbg-index] wrote:\n  {outp}.tdbg.meta\n  {outp}.tdbg.unitigs.tsv\n  {outp}.tdbg.k2u.tsv\n  {outp}.tids.txlen.tsv")

if __name__ == "__main__":
    main()
