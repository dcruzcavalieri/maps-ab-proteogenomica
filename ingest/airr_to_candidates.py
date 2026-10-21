"""
airr_to_candidates.py  —  AIRR (TRUST4 / pipeline do Gabriel) -> FASTA de candidatos
Lê um TSV AIRR e escreve um FASTA de região variável (VH/VL) para alimentar
run_integration.py. Usa sequence_aa / sequence_alignment_aa quando houver; senão
traduz sequence/sequence_vdj no frame indicado por v_sequence_start.

Ex.: python airr_to_candidates.py resultado_trust4_airr.tsv -o candidatos.fasta --min-consensus 2
"""
import argparse, csv, sys

CODON = {}
_b = "TCAG"; _a = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_i = 0
for a in _b:
    for b in _b:
        for c in _b:
            CODON[a+b+c] = _a[_i]; _i += 1

def translate(nt):
    nt = (nt or "").upper().replace("-", "").replace(".", "")
    p = "".join(CODON.get(nt[i:i+3], "X") for i in range(0, len(nt)-2, 3))
    return p.split("*")[0]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("airr"); ap.add_argument("-o", "--out", default="candidatos.fasta")
    ap.add_argument("--min-consensus", type=int, default=1)
    args = ap.parse_args()
    n = 0
    with open(args.airr) as fh, open(args.out, "w") as out:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            if r.get("productive", "T").upper().startswith("F"): continue
            try:
                if int(float(r.get("consensus_count", r.get("duplicate_count", 1)) or 1)) < args.min_consensus:
                    continue
            except ValueError:
                pass
            aa = (r.get("sequence_aa") or r.get("sequence_alignment_aa") or "").replace("-", "").replace(".", "")
            if not aa:
                aa = translate(r.get("sequence_vdj") or r.get("sequence") or "")
            if len(aa) < 60: continue
            sid = r.get("sequence_id", f"cand{n}")
            loc = r.get("locus", "?"); v = r.get("v_call", "?").split(",")[0]; j = r.get("j_call", "?").split(",")[0]
            out.write(f">{sid}|{loc}|{v}|{j}\n{aa}\n"); n += 1
    print(f"{n} candidatos -> {args.out}")

if __name__ == "__main__":
    main()
