"""
prep_ms_peptides.py  —  Normaliza tabelas de MS (PRIDE _psms.txt, MaxQuant,
MSFragger, PEAKS) para o formato de entrada de run_integration.py:
colunas  Peptide, ALC (%), Area.

Detecta automaticamente a coluna de sequência e remove modificações. Opcional:
--enzyme para renomear a saída (a protease é inferida do nome do arquivo em
run_integration.py, então inclua-a no nome: ..._Trypsin.csv, ..._Chymotrypsin.csv).

Ex.: python prep_ms_peptides.py 7649_JL_5a_psms.txt -o A1_S1_Elu_Trypsin.csv
"""
import argparse, os, re
import pandas as pd

SEQ_COLS = ["peptide", "sequence", "base sequence", "modified peptide",
            "peptide sequence", "annotated sequence", "pep"]
SCORE_COLS = ["alc (%)", "score", "ascore", "pep score", "andromeda score", "hyperscore", "-10lgp"]
AREA_COLS = ["area", "intensity", "precursor intensity", "ms1 area"]

def pick(cols, cands):
    low = {c.lower(): c for c in cols}
    for k in cands:
        if k in low: return low[k]
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile"); ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()
    df = pd.read_csv(args.infile, sep=None, engine="python")
    sc = pick(df.columns, SEQ_COLS)
    if not sc: raise SystemExit(f"coluna de peptídeo não encontrada em {args.infile}: {list(df.columns)[:12]}")
    pep = df[sc].astype(str).map(lambda s: re.sub(r"\([^)]*\)|\[[^\]]*\]|[^A-Za-z]", "", s).upper())
    out = pd.DataFrame({"Peptide": pep})
    scorec = pick(df.columns, SCORE_COLS)
    out["ALC (%)"] = pd.to_numeric(df[scorec], errors="coerce").fillna(90) if scorec else 90
    areac = pick(df.columns, AREA_COLS)
    out["Area"] = pd.to_numeric(df[areac], errors="coerce").fillna(1.0) if areac else 1.0
    out = out[out["Peptide"].str.len() >= 6].drop_duplicates()
    dst = args.out or os.path.splitext(args.infile)[0] + "_prep.csv"
    out.to_csv(dst, index=False)
    print(f"{len(out)} peptídeos -> {dst}  (col seq='{sc}', score='{scorec}', area='{areac}')")

if __name__ == "__main__":
    main()
