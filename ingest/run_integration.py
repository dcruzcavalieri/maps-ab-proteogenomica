"""
run_integration.py  —  Integração proteogenômica em DADOS REAIS
================================================================
Generaliza a POC (poc_insilico/poc_integracao.py) para entradas reais:

  * --candidates cands.fasta   candidatos de anticorpo do BCR-seq/RNA-seq
                               (VH/VL; p.ex. saída do TRUST4 / pipeline do Gabriel,
                               via ingest/airr_to_candidates.py).
  * --ms-dir DIR  |  --ms-csv F [F ...]
                               peptídeos de MS por protease/amostra. A protease é
                               inferida do nome do arquivo (Trypsin/Chymotrypsin/
                               LysC/AspN/GluC). Coluna aceita: Peptide|Sequence.

Saída (duas direções da validação mútua):
  (1) quais candidatos estão EXPRESSOS (peptídeos discriminantes + CDR3 corroborada);
  (2) CDR3 da MS ANCORADA no banco de candidatos (corrige erro de de novo).

Reusa o motor do MAPS-ab: maps_ab.cdr3_repertoire (FR4-anchored + cross-enzima).

Ex.:  python run_integration.py --candidates dados/candidatos_A1.fasta \
          --ms-dir dados/ms_peptides/A1  --out docs/integ_A1.json
"""
from __future__ import annotations
import argparse, glob, json, os, re, sys
from collections import defaultdict

# Resolve o pacote `maps_ab` (fork MAPS-ab-adaptado): env var MAPS_AB, senão
# pasta irmã ../../MAPS-ab-adaptado. Instale o fork ou defina MAPS_AB.
_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (os.environ.get("MAPS_AB"), os.path.join(_HERE, "..", "..", "MAPS-ab-adaptado")):
    if _p and os.path.isdir(os.path.join(_p, "maps_ab")):
        sys.path.insert(0, os.path.abspath(_p)); break
try:
    from maps_ab.cdr3_repertoire import extract_records, build_repertoire, enzyme_of, _canon
except ImportError:
    sys.exit("ERRO: pacote 'maps_ab' não encontrado. Defina MAPS_AB=/caminho/para/MAPS-ab-adaptado "
             "ou coloque o fork como pasta irmã deste repositório.")
import pandas as pd

ENZYMES = ["Trypsin", "Chymotrypsin", "LysC", "AspN", "GluC"]
_FR4 = re.compile(r"[FY][YFH]C(.{4,30}?)([WF])G[A-Z]G[TS]")   # CDR3 de uma cadeia


def cleave(seq, enzyme, max_missed=2):
    if enzyme == "Trypsin":        s = [i for i in range(len(seq)-1) if seq[i] in "KR" and seq[i+1] != "P"]
    elif enzyme == "Chymotrypsin": s = [i for i in range(len(seq)-1) if seq[i] in "FYW" and seq[i+1] != "P"]
    elif enzyme == "LysC":         s = [i for i in range(len(seq)-1) if seq[i] == "K"]
    elif enzyme == "AspN":         s = [i-1 for i in range(1, len(seq)) if seq[i] == "D"]
    elif enzyme == "GluC":         s = [i for i in range(len(seq)-1) if seq[i] in "ED"]
    else: s = []
    cuts = [-1] + sorted(set(x for x in s if 0 <= x < len(seq))) + [len(seq)-1]
    out = []
    for i in range(len(cuts)-1):
        for m in range(max_missed+1):
            if i+1+m < len(cuts):
                out.append(seq[cuts[i]+1: cuts[i+1+m]+1])
    return [p for p in out if 7 <= len(p) <= 40]


def candidate_cdr3(seq):
    m = _FR4.search(seq)
    return m.group(1) if m else None


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs


def load_ms(files):
    """-> dict enzyme -> DataFrame(Peptide, ALC (%), Area) ; lista (pep,enz)."""
    per_enz = {}
    for f in files:
        try:
            df = pd.read_csv(f, sep=None, engine="python")
        except Exception as e:
            print(f"  (skip {os.path.basename(f)}: {e})"); continue
        col = next((c for c in df.columns if c.lower() in
                    ("peptide", "sequence", "peptide sequence", "pep")), None)
        if not col:
            print(f"  (sem coluna de peptídeo em {os.path.basename(f)})"); continue
        out = pd.DataFrame({"Peptide": df[col].astype(str)})
        out["ALC (%)"] = pd.to_numeric(df.get("ALC (%)", df.get("Score", 90)), errors="coerce").fillna(90)
        out["Area"] = pd.to_numeric(df.get("Area", df.get("Intensity", 1.0)), errors="coerce").fillna(1.0)
        enz = enzyme_of(f)
        per_enz.setdefault(enz, []).append(out)
    return {e: pd.concat(v, ignore_index=True) for e, v in per_enz.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True, help="FASTA de candidatos VH/VL (banco)")
    ap.add_argument("--ms-dir", default=None)
    ap.add_argument("--ms-csv", nargs="*", default=[])
    ap.add_argument("--min-alc", type=float, default=70.0)
    ap.add_argument("--min-discrim", type=int, default=3, help="peptídeos discriminantes p/ chamar expresso")
    ap.add_argument("--out", default="integracao_resultado.json")
    args = ap.parse_args()

    cands = read_fasta(args.candidates)
    files = list(args.ms_csv) + (sorted(glob.glob(os.path.join(args.ms_dir, "**", "*.csv"),
                                                   recursive=True)) +
                                 sorted(glob.glob(os.path.join(args.ms_dir, "**", "*.txt"),
                                                   recursive=True)) if args.ms_dir else [])
    if not files: raise SystemExit("Sem CSVs de MS (use --ms-dir ou --ms-csv).")
    ms = load_ms(files)
    enz_present = sorted(ms.keys())

    # peptídeos teóricos (canon) por candidato
    cand_pep = {cid: set(_canon(p) for e in ENZYMES for p in cleave(seq, e))
                for cid, seq in cands.items()}
    cand_cdr3 = {cid: candidate_cdr3(seq) for cid, seq in cands.items()}

    # peptídeos observados (ALC>=min)
    obs = []
    for e, df in ms.items():
        for p, a in zip(df["Peptide"], df["ALC (%)"]):
            if a >= args.min_alc:
                obs.append((re.sub(r"\([^)]*\)", "", str(p)).strip(), e))

    def approx(o, pepset, mm=1):
        o = _canon(o); L = len(o)
        for q in pepset:
            if len(q) == L and sum(1 for x, y in zip(o, q) if x != y) <= mm:
                return True
        return False

    support = defaultdict(lambda: {"discrim": 0, "shared": 0, "enz": set()})
    for o, e in obs:
        hits = [cid for cid in cands if approx(o, cand_pep[cid])]
        if len(hits) == 1:
            support[hits[0]]["discrim"] += 1; support[hits[0]]["enz"].add(e)
        elif len(hits) > 1:
            for cid in hits: support[cid]["shared"] += 1

    # CDR3 recuperadas da MS pelo motor do MAPS-ab
    recs = []
    for e, df in ms.items():
        recs += extract_records(df, args.min_alc, 5, 26, enzyme=e, source=e)
    clones = build_repertoire(recs, min_enzymes=1)

    def cdr3_confirmed(cid):
        ref = cand_cdr3[cid]
        if not ref: return (False, 0)
        rc = _canon(ref); best = (0.0, 0)
        for c in clones:
            cc = _canon(c.cdr3); n = min(len(cc), len(rc))
            if n < 6: continue
            idv = sum(x == y for x, y in zip(cc[-n:], rc[-n:]))/n
            if idv > best[0]: best = (idv, c.n_enzymes)
        return (best[0] >= 0.8, best[1])

    rows = []
    for cid in cands:
        disc = support[cid]["discrim"]; nenz = len(support[cid]["enz"])
        conf, cne = cdr3_confirmed(cid)
        rows.append({"candidate": cid, "discrim_peptides": disc, "enzymes": nenz,
                     "cdr3_confirmed": conf, "cdr3_n_enzymes": cne,
                     "expressed_call": (disc >= args.min_discrim) or conf,
                     "cdr3": cand_cdr3[cid]})
    rows.sort(key=lambda r: (-int(r["expressed_call"]), -r["discrim_peptides"]))

    # (2) ancoramento das CDR3 da MS no banco
    db_cdr3 = {cid: _canon(c) for cid, c in cand_cdr3.items() if c}
    def cterm(a, b):
        n = min(len(a), len(b)); return sum(x == y for x, y in zip(a[-n:], b[-n:]))/n if n else 0
    anchored = []
    for c in clones:
        raw = _canon(c.cdr3)
        if db_cdr3:
            best_cid = max(db_cdr3, key=lambda k: cterm(raw, db_cdr3[k]))
            anchored.append({"ms_cdr3": c.cdr3, "n_enzymes": c.n_enzymes,
                             "anchored_to": best_cid, "anchored_cdr3": cand_cdr3[best_cid],
                             "cterm_id": round(cterm(raw, db_cdr3[best_cid]), 3)})

    report = {"candidates_n": len(cands), "ms_files": len(files), "enzymes_present": enz_present,
              "observed_peptides": len(obs),
              "expressed_calls": sum(1 for r in rows if r["expressed_call"]),
              "per_candidate": rows, "anchored_cdr3": anchored[:50]}
    json.dump(report, open(args.out, "w"), indent=2, ensure_ascii=False)

    print(f"Candidatos: {len(cands)} | MS: {len(obs)} peptídeos | enzimas: {enz_present}")
    print(f"{'candidato':28}{'disc':>5}{'enz':>4}{'CDR3conf':>9}  chamada")
    for r in rows[:30]:
        print(f"{r['candidate'][:28]:28}{r['discrim_peptides']:>5}{r['enzymes']:>4}"
              f"{('sim('+str(r['cdr3_n_enzymes'])+'e)') if r['cdr3_confirmed'] else '-':>9}  "
              f"{'EXPRESSO' if r['expressed_call'] else '-'}")
    print(f"\nExpressos: {report['expressed_calls']}/{len(cands)}  ->  {args.out}")


if __name__ == "__main__":
    main()
