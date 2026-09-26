import glob, os, random, sys, json
sys.path.insert(0, "/home/claude/scss")
from stp_parse import parse_stp
from directionize2 import directionize_q
OUT = "/home/claude/scss/directed2"
Q = {"A": 1.0, "B": 0.0, "C": 0.5}
SEEDS = [42, 1, 2]
files = sorted(glob.glob("/home/claude/scss/scss-benchmark-instances/*/*.stp"))
meta = []
for f in files:
    inst = parse_stp(f); setn = f.split("/")[-2]; base = inst.name.replace(".stp", ".txt")
    for lab, q in Q.items():
        for seed in (SEEDS if lab != "A" else [42]):
            d, nbq, nbfix, now = directionize_q(inst, q, random.Random(seed))
            d = sorted(set(d))
            dd = os.path.join(OUT, f"{lab}_s{seed}", setn); os.makedirs(dd, exist_ok=True)
            with open(os.path.join(dd, base), "w") as w:
                w.write(f"{inst.n} {len(d)} {len(inst.terminals)}\n")
                for u, v in d: w.write(f"{u} {v}\n")
                w.write(" ".join(map(str, inst.terminals)) + "\n")
            meta.append(dict(set=setn, name=base, lab=lab, seed=seed, m_und=len(inst.edges),
                             bidir_q=nbq, bidir_fix=nbfix, oneway=now))
json.dump(meta, open(OUT + "/meta.json", "w"), indent=1)
print("exported", len(meta))
