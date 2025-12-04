
from python import sys
from common import read_fastq, effective_length, write_abundance

class EqClass:
    def __init__(self, members, count):
        self.members = members
        self.count = float(count)

def load_tdbg(prefix: str):
    meta_path    = prefix + ".tdbg.meta"
    k2u_path     = prefix + ".tdbg.k2u.tsv"
    unitig_path  = prefix + ".tdbg.unitigs.tsv"
    tx_path      = prefix + ".tids.txlen.tsv"

    k = 0
    with open(meta_path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("k="):
                left = s.split("\t")[0]
                k = int(left.split("=")[1])
                break

    tids, lengths, tid2ix = [], [], {}
    with open(tx_path) as f:
        for line in f:
            s = line.strip()
            if not s or "\t" not in s: continue
            tid, Ls = s.split("\t")
            tid2ix[tid] = len(tids)
            tids.append(tid)
            lengths.append(int(Ls))

    unitig2tix = {}
    with open(unitig_path) as f:
        for line in f:
            s = line.strip()
            if not s or "\t" not in s: continue
            uid_str, tids_str = s.split("\t")
            uid = int(uid_str)
            tx_ids = [t for t in tids_str.split(",") if t != ""]
            seen, lst = {}, []
            for t in tx_ids:
                if t in tid2ix:
                    ix = tid2ix[t]
                    if ix not in seen:
                        seen[ix] = 1
                        lst.append(ix)
            lst.sort()
            unitig2tix[uid] = lst

    k2tix = {}
    with open(k2u_path) as f:
        for line in f:
            s = line.strip()
            if not s or "\t" not in s: continue
            km, uid_str = s.split("\t")
            uid = int(uid_str)
            lst = unitig2tix.get(uid, [])
            if lst:
                k2tix[km] = lst

    return k, k2tix, tids, lengths

# def pseudoalign_read_on_tdbg(seq: str, k: int, k2tix: dict[str, list[int]],
#                              mode: str, min_hits: int, skip_budget: int) -> list[int]:
#     n = len(seq)
#     hits = []
#     for i in range(n - k + 1):
#         km = seq[i:i+k]
#         if km in k2tix:
#             hits.append(k2tix[km])
#     if len(hits) < min_hits:
#         return []

#     if mode == "union":
#         mark, merged = {}, []
#         for h in hits:
#             for ix in h:
#                 if ix not in mark:
#                     mark[ix] = 1
#                     merged.append(ix)
#         merged.sort()
#         return merged
#     else:
#         s = hits[0][:]
#         s.sort()

#         def _inter(a, b):
#             a2, b2 = a[:], b[:]
#             a2.sort(); b2.sort()
#             i = 0; j = 0; out = []
#             while i < len(a2) and j < len(b2):
#                 if a2[i] == b2[j]:
#                     out.append(a2[i]); i += 1; j += 1
#                 elif a2[i] < b2[j]:
#                     i += 1
#                 else:
#                     j += 1
#             return out

#         for h in hits[1:]:
#             s = _inter(s, h)
#             if not s: break
#         return s

def pseudoalign_read_on_tdbg(seq: str, k: int, k2tix: dict[str, list[int]],
                             mode: str, min_hits: int, skip_budget: int) -> list[int]:
    """Pseudoalign one read using a simple skip-budget heuristic.

    skip_budget = how many consecutive non-matching k-mers we are allowed to *check*
    before we start skipping whole k-mer blocks.
    """
    n = len(seq)
    hits: list[list[int]] = []

    # اگر skip_budget <= 0 باشه، مثل نسخه‌ی قدیمی رفتار می‌کنیم (هیچ اسکیپی نداریم)
    if skip_budget is None or skip_budget <= 0:
        i = 0
        while i <= n - k:
            km = seq[i:i+k]
            if km in k2tix:
                hits.append(k2tix[km])
            i += 1
    else:
        i = 0
        remaining = skip_budget
        while i <= n - k:
            km = seq[i:i+k]
            if km in k2tix:
                hits.append(k2tix[km])
                remaining = skip_budget
                i += 1
            else:
                if remaining > 0:
                    remaining -= 1
                    i += 1
                else:
                    i += k
                    remaining = skip_budget

    if len(hits) < min_hits:
        return []

    if mode == "union":
        mark, merged = {}, []
        for h in hits:
            for ix in h:
                if ix not in mark:
                    mark[ix] = 1
                    merged.append(ix)
        merged.sort()
        return merged
    else:
        s = hits[0][:]
        s.sort()

        def _inter(a, b):
            a2, b2 = a[:], b[:]
            a2.sort(); b2.sort()
            i = 0; j = 0; out = []
            while i < len(a2) and j < len(b2):
                if a2[i] == b2[j]:
                    out.append(a2[i]); i += 1; j += 1
                elif a2[i] < b2[j]:
                    i += 1
                else:
                    j += 1
            return out

        for h in hits[1:]:
            s = _inter(s, h)
            if not s:
                break
        return s

def build_equivalence_classes(reads1: list[str], reads2,
                              k: int, k2tix: dict[str, list[int]],
                              mode: str, min_hits: int, skip_budget: int) -> dict[str, EqClass]:
    ecs, n = {}, len(reads1)
    paired = reads2 is not None

    for i in range(n):
        s1 = pseudoalign_read_on_tdbg(reads1[i], k, k2tix, mode, min_hits, skip_budget)
        if not s1: continue
        comb = s1

        if paired and i < len(reads2):
            s2 = pseudoalign_read_on_tdbg(reads2[i], k, k2tix, mode, min_hits, skip_budget)
            if mode == "union":
                if s2:
                    seen, merged = {}, []
                    for x in comb:
                        if x not in seen: seen[x] = 1; merged.append(x)
                    for x in s2:
                        if x not in seen: seen[x] = 1; merged.append(x)
                    merged.sort()
                    comb = merged
            else:
                if not s2: continue
                a, b = comb[:], s2[:]
                a.sort(); b.sort()
                i1 = 0; j1 = 0; inter = []
                while i1 < len(a) and j1 < len(b):
                    if a[i1] == b[j1]:
                        inter.append(a[i1]); i1 += 1; j1 += 1
                    elif a[i1] < b[j1]:
                        i1 += 1
                    else:
                        j1 += 1
                if not inter: continue
                comb = inter

        if not comb: continue
        key = ",".join(str(x) for x in comb)
        if key in ecs: ecs[key].count += 1.0
        else:          ecs[key] = EqClass(comb, 1.0)
    return ecs

def em_quant(classes: dict[str, EqClass], n_tx: int, lengths: list[int],
             frag_len: int, max_iter: int = 200, tol: float = 1e-6):
    eff = [effective_length(lengths[i], frag_len) for i in range(n_tx)]
    Z = sum(eff);  Z = Z if Z > 0.0 else 1.0
    theta = [eff[i] / Z for i in range(n_tx)]

    for _ in range(max_iter):
        alloc = [0.0] * n_tx
        for ec in classes.values():
            denom = 0.0
            for ix in ec.members:
                denom += theta[ix] * eff[ix]
            if denom <= 0.0: continue
            invd = 1.0 / denom
            for ix in ec.members:
                alloc[ix] += ec.count * (theta[ix] * eff[ix]) * invd

        total = sum(alloc)
        new_theta = [a / total for a in alloc] if total > 0.0 else theta[:]
        diff = 0.0
        for i in range(n_tx):
            d = new_theta[i] - theta[i]
            if d < 0.0: d = -d
            diff += d
        theta = new_theta
        if diff < tol: break

    est_counts = [0.0] * n_tx
    for ec in classes.values():
        denom = 0.0
        for ix in ec.members:
            denom += theta[ix] * eff[ix]
        if denom <= 0.0: continue
        invd = 1.0 / denom
        for ix in ec.members:
            est_counts[ix] += ec.count * (theta[ix] * eff[ix]) * invd

    rates = [ (est_counts[i] / eff[i]) if eff[i] > 0.0 else 0.0 for i in range(n_tx) ]
    sum_rates = sum(rates);  scale = (sum_rates / 1_000_000.0) if sum_rates > 0.0 else 1.0
    tpm = [ (r / scale) if scale > 0.0 else 0.0 for r in rates ]
    return est_counts, tpm, eff

_rng_state = 88172645463393265
def rng_seed(seed: int):
    global _rng_state
    if seed <= 0: seed = 88172645463393265
    _rng_state = seed & ((1<<64)-1)
    if _rng_state == 0: _rng_state = 88172645463393265

def rng_next_u64() -> int:
    global _rng_state
    x = _rng_state
    x ^= (x >> 12) & ((1<<64)-1)
    x ^= (x << 25) & ((1<<64)-1)
    x ^= (x >> 27) & ((1<<64)-1)
    _rng_state = x
    return (x * 2685821657736338717) & ((1<<64)-1)

def rng_u01() -> float:
    return (rng_next_u64() >> 11) * (1.0 / (1<<53))

def pois(lam: float) -> int:
    if lam <= 0.0: return 0
    if lam > 50.0:
        import math
        u1 = rng_u01(); u2 = rng_u01()
        z = ((-2.0 * math.log(max(u1,1e-12))) ** 0.5) * (math.cos(2.0*3.141592653589793*u2))
        val = int(lam + z * (lam ** 0.5))
        return val if val > 0 else 0
    import math
    L = math.exp(-lam)
    k = 0; p = 1.0
    while True:
        k += 1
        u = rng_u01()
        p *= u
        if p <= L: break
    return k - 1

def bootstrap_once(classes: dict[str, EqClass], n_tx: int, lengths: list[int], frag_len: int) -> list[float]:
    b_classes = {}
    for key, ec in classes.items():
        c = pois(ec.count)
        if c > 0:
            b_classes[key] = EqClass(ec.members, float(c))
    est_counts, tpm, _ = em_quant(b_classes, n_tx, lengths, frag_len)
    return tpm

def summarize_bootstrap(tpms: list[list[float]]):
    import math
    nB = len(tpms)
    if nB == 0: return []
    n_tx = len(tpms[0])
    means = [0.0]*n_tx
    for b in range(nB):
        row = tpms[b]
        for i in range(n_tx):
            means[i] += row[i]
    for i in range(n_tx): means[i] /= float(nB)

    sds = [0.0]*n_tx
    for b in range(nB):
        row = tpms[b]
        for i in range(n_tx):
            d = row[i] - means[i]
            sds[i] += d*d
    for i in range(n_tx): sds[i] = (sds[i] / float(nB - 1)) ** 0.5 if nB > 1 else 0.0

    ses = [ (sds[i] / (float(nB) ** 0.5)) if nB > 0 else 0.0 for i in range(n_tx) ]
    cvs = [ (sds[i]/means[i]) if means[i] > 0.0 else 0.0 for i in range(n_tx) ]
    ciL = [ means[i] - 1.96 * ses[i] for i in range(n_tx) ]
    ciH = [ means[i] + 1.96 * ses[i] for i in range(n_tx) ]
    return [(means[i], sds[i], ses[i], cvs[i], ciL[i], ciH[i], nB) for i in range(n_tx)]

def write_bootstrap_summary(path: str, tids: list[str], stats_list):
    with open(path, "w") as w:
        w.write("target_id\tmean_tpm\tsd_tpm\tse_tpm\tcv_tpm\tci95_low\tci95_high\tB_eff\n")
        for i in range(len(tids)):
            m, sd, se, cv, lo, hi, B = stats_list[i]
            w.write(f"{tids[i]}\t{m:.6f}\t{sd:.6f}\t{se:.6f}\t{cv:.6f}\t{lo:.6f}\t{hi:.6f}\t{B}\n")

def _count_reads_with_hits(reads: list[str], k: int, k2tix: dict[str, list[int]]) -> int:
    c = 0
    for s in reads:
        n = len(s); ok = False
        for j in range(n - k + 1):
            km = s[j:j+k]
            if km in k2tix and len(k2tix[km]) > 0:
                ok = True; break
        if ok: c += 1
    return c

def _cli_args() -> list[str]:
    out = []
    try:
        if hasattr(sys, "argv") and sys.argv and len(sys.argv) > 1:
            for a in sys.argv[1:]:
                s = f"{a}"
                if s and s != "--": out.append(s)
    except: pass
    try:
        if hasattr(sys, "args") and sys.args and len(sys.args) > 0:
            arr = [f"{a}" for a in sys.args]
            if arr and (arr[0].endswith(".py") or arr[0].endswith(".codon")):
                arr = arr[1:]
            for s in arr:
                if s and s != "--": out.append(s)
    except: pass
    try:
        from python import os as pyos
        envs = pyos.getenv("CODON_ARGS")
        if envs is not None:
            for tok in f"{envs}".strip().split():
                if tok != "--": out.append(f"{tok}")
    except: pass
    return out

def main():
    args = _cli_args()
    if len(args) < 4:
        print("Usage:\n  codon run quant_dbg.py -- <tdbg_prefix> <frag_len> <reads_1.fastq> [reads_2.fastq] <out_dir> [--mode union|intersect] [--min-khits N] [--skip-budget S] [--bootstrap B] [--seed S]")
        print("Args seen:", args)
        sys.exit(1)

    idx_prefix = args[0].strip()
    frag_len   = int(args[1])
    reads1_path = args[2].strip()

    i = 3
    reads2_path = None
    if i < len(args) and not args[i].startswith("--"):
        reads2_path = args[i].strip(); i += 1
    if i >= len(args):
        print("Missing output directory!")
        sys.exit(1)
    out_dir = args[i].strip(); i += 1  # (kept for CLI) we write in CWD

    mode = "union"
    min_hits = 1
    skip_budget = 3
    B = 0
    seed = 123
    while i < len(args):
        a = args[i]
        if a == "--mode" and i + 1 < len(args):
            mode = args[i + 1]; i += 2
        elif a == "--min-khits" and i + 1 < len(args):
            min_hits = int(args[i + 1]); i += 2
        elif a == "--skip-budget" and i + 1 < len(args):
            skip_budget = int(args[i + 1]); i += 2
        elif a == "--bootstrap" and i + 1 < len(args):
            B = int(args[i + 1]); i += 2
        elif a == "--seed" and i + 1 < len(args):
            seed = int(args[i + 1]); i += 2
        else:
            i += 1

    k, k2tix, tids, lengths = load_tdbg(idx_prefix)
    reads1 = read_fastq(reads1_path)
    reads2 = read_fastq(reads2_path) if reads2_path else None

    print("[debug] SE1 reads with >=1 hit:", _count_reads_with_hits(reads1, k, k2tix))
    if reads2 is not None:
        print("[debug] SE2 reads with >=1 hit:", _count_reads_with_hits(reads2, k, k2tix))

    ecs = build_equivalence_classes(reads1, reads2, k, k2tix, mode, min_hits, skip_budget)
    print(f"[quant-dbg] equivalence classes = {len(ecs)}")

    est_counts, tpm, eff = em_quant(ecs, len(tids), lengths, frag_len)

    outp = "./abundance.tsv"
    rows = [(tids[i], lengths[i], eff[i], est_counts[i], tpm[i]) for i in range(len(tids))]
    print(f"[debug] writing to <{outp}>")
    try:
        write_abundance(outp, rows)
        print(f"[quant-dbg] wrote → {outp}")
    except Exception as e:
        outp2 = "./abundance_out.tsv"
        print(f"[warn] failed to write <{outp}>; fallback → <{outp2}>; err={e}")
        write_abundance(outp2, rows)
        print(f"[quant-dbg] wrote → {outp2}")

    if B > 0:
        rng_seed(seed)
        tpm_mat = []
        for _ in range(B):
            tpm_b = bootstrap_once(ecs, len(tids), lengths, frag_len)
            tpm_mat.append(tpm_b)
        stats_list = summarize_bootstrap(tpm_mat)
        bs_path = "./abundance_bootstrap_summary.tsv"
        try:
            write_bootstrap_summary(bs_path, tids, stats_list)
            print(f"[quant-dbg] bootstrap summary → {bs_path}")
        except Exception as e:
            bs2 = "./abundance_bootstrap_summary_out.tsv"
            print(f"[warn] failed to write <{bs_path}>; fallback → <{bs2}>; err={e}")
            write_bootstrap_summary(bs2, tids, stats_list)
            print(f"[quant-dbg] bootstrap summary → {bs2}")

if __name__ == "__main__":
    main()
