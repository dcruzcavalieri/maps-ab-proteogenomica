"""
poc_integracao.py  —  Prova de conceito IN-SILICO da integração proteogenômica
==============================================================================
Demonstra, com RESPOSTA CONHECIDA, a ideia central do braço MAPS-ab x TCC:

  RNA-seq (candidatos de sequência)  +  MS de novo (o que está expresso)
  → validação mútua.

Modelo:
  * BANCO (candidatos de RNA-seq): os 5 candidatos reconstruídos pelo pipeline do
    Gabriel + F59 (HC/LC) + distratores (anticorpos humanos reais, não relacionados).
  * SORO (verdade-terreno): só um subconjunto está EXPRESSO como proteína.
  * MS simulada: digestão multi-protease (tripsina/quimotripsina/AspN) das cadeias
    EXPRESSAS, com cobertura parcial, ambiguidade I/L, erros de de novo e clivagem
    perdida — gerando "leituras de novo" como as do PEAKS que o MAPS-ab consome.
  * INTEGRAÇÃO:
      (1) MS confirma quais candidatos do banco estão expressos (peptídeos
          discriminantes + CDR3 corroborada, via motor do MAPS-ab);
      (2) o banco ANCORA/corrige as CDR3 ruidosas da MS (de novo puro vs ancorado).

Robustez: a simulação é repetida em N sementes; reporta média (min–max).
Saídas: poc_resultado.json, ms_sim/*.csv (semente de exemplo).
"""
from __future__ import annotations
import json, os, random, sys
from collections import defaultdict

# Resolve o pacote `maps_ab` (fork MAPS-ab-adaptado): variável de ambiente MAPS_AB,
# senão pasta irmã ../../MAPS-ab-adaptado. Instale o fork ou aponte MAPS_AB p/ ele.
_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (os.environ.get("MAPS_AB"), os.path.join(_HERE, "..", "..", "MAPS-ab-adaptado")):
    if _p and os.path.isdir(os.path.join(_p, "maps_ab")):
        sys.path.insert(0, os.path.abspath(_p)); break
try:
    from maps_ab.cdr3_repertoire import extract_records, build_repertoire, _canon
except ImportError:
    sys.exit("ERRO: pacote 'maps_ab' não encontrado. Defina MAPS_AB=/caminho/para/MAPS-ab-adaptado "
             "ou coloque o fork como pasta irmã deste repositório.")
import pandas as pd

OUT = _HERE
ENZ = ["Trypsin", "Chymotrypsin", "AspN"]
P_OBS = 0.60      # cobertura parcial (fração de peptídeos teóricos observados)
P_SUB = 0.08      # taxa de erro de de novo (miscall por resíduo) — realista
N_SEEDS = 12
AA = "ACDEFGHIKLMNPQRSTVWY"

# -----------------------------------------------------------------------------
# BANCO de candidatos (RNA-seq) + F59 + distratores  (fixo)
# -----------------------------------------------------------------------------
def load_DB():
    cands = json.load(open(os.path.join(_HERE, "candidatos_pipeline.json")))
    DB = {}
    for c in cands:
        DB[c["id"]] = {"chain": c["chain"], "aa": c["aa"],
                       "cdr3": c["cdr3"].strip("C"), "source": "RNAseq_cand"}
    DB["F59_HC"] = {"chain":"heavy","aa":"EPELVESGGLAQPGTSLRLSCEAASGFTFDDYAMHWVRQAPGRALEWVSGISWSSDNIAYSDSVEGRFTISRDNAKNSLYLQMNSLRLDDTAFYYCAKDVRPYYDFWAFDSWGRGTPVTV","cdr3":"AKDVRPYYDFWAFDS","source":"F59_ref"}
    DB["F59_LC"] = {"chain":"light","aa":"QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYYVSWYQQHPGKAPKLLISEVSDRPSGVSRFSGSKSGNTASLTISGLQAEDESMYFCSSYTDLTSFVVFGGGTKLTV","cdr3":"SSYTDLTSFVV","source":"F59_ref"}
    DISTR = {
     "2D22_HC":("heavy","EVQLVESGAEVKKPGSSVKVSCKASGGTFNNYAISWVRQAPGQGLEWMGGIIPIFGGANYAQKFQGRVTITADRSTSTVYMELSGLRSEDTAVYYCARRPQSIFDWNFDLWGRGTLVTVSSAGTKGPS"),
     "2D22_LC":("light","QSVLTQPPSASGTPGQRVTISCSGSSSNVGSNYVYWYQQLPGTAPKLLIYRNNRRPSGVPDRFSGSKSGTSASLAISGLRSEDEADYYCATWDDSLSGLVFGGGTKLTVLGQPKA"),
     "EDE1C8_HC":("heavy","EVQLVESGGGLVQPGGSLRLSCSASGFTFSTYSMHWVRQAPGKGLEYVSAITGEGDSAFYADSVKGRFTISRDNSKNTLYFEMNSLRPEDTAVYYCVGGYSNFYYYYTMDVWGQGTTVTVSS"),
     "EDE1C8_LC":("light","RSEIVLTQSPATLSLSPGERATLSCRASQSISTFLAWYQHKPGQAPRLLIYDASTRATGVPARFSGSRSGTDFTLTISTLEPEDFAVYYCQQRYNWPPYTFGQGTKVEIK"),
    }
    for k,(ch,aa) in DISTR.items():
        DB[k] = {"chain":ch, "aa":aa, "cdr3":None, "source":"distrator"}
    return DB

DB = load_DB()
EXPRESSED = {"assemble181_0", "assemble95_0", "F59_HC", "F59_LC"}   # verdade-terreno

# -----------------------------------------------------------------------------
# digestão in-silico e matching
# -----------------------------------------------------------------------------
def cleave(seq, enzyme, max_missed=2):
    if enzyme=="Trypsin":        sites=[i for i in range(len(seq)-1) if seq[i] in "KR" and seq[i+1]!="P"]
    elif enzyme=="Chymotrypsin": sites=[i for i in range(len(seq)-1) if seq[i] in "FYW" and seq[i+1]!="P"]
    elif enzyme=="AspN":         sites=[i-1 for i in range(1,len(seq)) if seq[i]=="D"]
    else: sites=[]
    cuts=[-1]+sorted(set(s for s in sites if 0<=s<len(seq)))+[len(seq)-1]
    peps=[]
    for i in range(len(cuts)-1):
        for m in range(0,max_missed+1):
            if i+1+m<len(cuts): peps.append(seq[cuts[i]+1:cuts[i+1+m]+1])
    return [p for p in peps if 7<=len(p)<=35]

def theoretical(seq):
    s=set()
    for e in ENZ:
        for p in cleave(seq,e): s.add(_canon(p))
    return s
DB_PEP = {cid: theoretical(v["aa"]) for cid,v in DB.items()}

def approx_in(obs, pepset, max_mm=1):
    o=_canon(obs); L=len(o)
    for q in pepset:
        if len(q)!=L: continue
        if sum(1 for x,y in zip(o,q) if x!=y)<=max_mm: return True
    return False

def cterm_id(a,b):
    n=min(len(a),len(b)); return sum(x==y for x,y in zip(a[-n:],b[-n:]))/n if n else 0.0

# -----------------------------------------------------------------------------
# uma rodada completa
# -----------------------------------------------------------------------------
def run_once(seed, write_csv=False):
    rnd=random.Random(seed)
    # MS simulada (só expressas)
    ms={e:[] for e in ENZ}
    for cid in EXPRESSED:
        for e in ENZ:
            for pep in set(cleave(DB[cid]["aa"],e)):
                if rnd.random()<P_OBS:
                    obs="".join(rnd.choice(AA) if rnd.random()<P_SUB else a for a in pep)
                    ms[e].append((obs, rnd.randint(72,99)))
    n_ms=sum(len(v) for v in ms.values())
    if write_csv:
        os.makedirs(os.path.join(OUT,"ms_sim"),exist_ok=True)
        for e in ENZ:
            pd.DataFrame({"Peptide":[p for p,_ in ms[e]],"ALC (%)":[a for _,a in ms[e]],
                          "Area":[rnd.randint(10**5,10**8) for _ in ms[e]]}
                         ).to_csv(os.path.join(OUT,"ms_sim",f"MSsim_{e}.csv"),index=False)

    # (1) atribuição de peptídeos ao BANCO
    support=defaultdict(lambda:{"discrim":0,"enz":set()})
    for e in ENZ:
        for obs,_ in ms[e]:
            hits=[cid for cid in DB if approx_in(obs,DB_PEP[cid])]
            if len(hits)==1:
                support[hits[0]]["discrim"]+=1; support[hits[0]]["enz"].add(e)
    # CDR3 recuperadas pela MS (motor MAPS-ab)
    recs=[]
    for e in ENZ:
        df=pd.DataFrame({"Peptide":[p for p,_ in ms[e]],"ALC (%)":[a for _,a in ms[e]],"Area":[1.0]*len(ms[e])})
        recs+=extract_records(df,70,5,26,enzyme=e,source=e)
    clones=build_repertoire(recs,min_enzymes=1)
    def cdr3_conf(cid):
        ref=DB[cid]["cdr3"]
        if not ref: return (False,0)
        rc=_canon(ref); best=(0,0)
        for c in clones:
            cc=_canon(c.cdr3); n=min(len(cc),len(rc))
            if n<6: continue
            idv=sum(x==y for x,y in zip(cc[-n:],rc[-n:]))/n
            if idv>best[0]: best=(idv,c.n_enzymes)
        return (best[0]>=0.8, best[1])
    rows=[]
    for cid,v in DB.items():
        disc=support[cid]["discrim"]; conf,ne=cdr3_conf(cid)
        rows.append({"id":cid,"source":v["source"],"truth":cid in EXPRESSED,
                     "discrim":disc,"enz":len(support[cid]["enz"]),
                     "cdr3_conf":conf,"cdr3_nenz":ne,"call":(disc>=3) or conf})
    tp=sum(1 for r in rows if r["call"] and r["truth"]); fp=sum(1 for r in rows if r["call"] and not r["truth"])
    fn=sum(1 for r in rows if not r["call"] and r["truth"]); tn=sum(1 for r in rows if not r["call"] and not r["truth"])
    prec=tp/(tp+fp) if tp+fp else 1.0; rec=tp/(tp+fn) if tp+fn else 1.0
    f1=2*prec*rec/(prec+rec) if prec+rec else 0.0

    # (2) ancoramento das CDR3 pelo BANCO
    true_cdr3={cid:_canon(DB[cid]["cdr3"]) for cid in EXPRESSED if DB[cid]["cdr3"]}
    db_cdr3={cid:_canon(v["cdr3"]) for cid,v in DB.items() if v["cdr3"]}
    ev=[]
    for c in clones:
        raw=_canon(c.cdr3)
        tgt=max(true_cdr3.values(),key=lambda t:cterm_id(raw,t))
        id_raw=cterm_id(raw,tgt)
        anc=max(db_cdr3.values(),key=lambda t:cterm_id(raw,t))
        id_anc=cterm_id(anc,tgt)
        if id_raw>=0.6 or id_anc>=0.6: ev.append((id_raw,id_anc))
    mr=sum(a for a,_ in ev)/len(ev) if ev else 0; ma=sum(b for _,b in ev)/len(ev) if ev else 0
    return {"n_ms":n_ms,"prec":prec,"rec":rec,"f1":f1,"tp":tp,"fp":fp,"fn":fn,"tn":tn,
            "mr":mr,"ma":ma,"rows":rows,"n_ev":len(ev)}

# -----------------------------------------------------------------------------
# múltiplas sementes
# -----------------------------------------------------------------------------
res=[run_once(s, write_csv=(s==0)) for s in range(N_SEEDS)]
def agg(key):
    xs=[r[key] for r in res]; return (sum(xs)/len(xs), min(xs), max(xs))
precs,recs_,f1s=agg("prec"),agg("rec"),agg("f1")
mr,ma=agg("mr"),agg("ma")
nms=agg("n_ms")
detailed=res[0]

report={
 "config":{"proteases":ENZ,"cobertura_obs":P_OBS,"erro_denovo":P_SUB,"n_sementes":N_SEEDS,
           "banco_n":len(DB),"expressas_verdade":sorted(EXPRESSED)},
 "classificacao_expressao_media":{"precisao":round(precs[0],3),"recall":round(recs_[0],3),
     "f1":round(f1s[0],3),"f1_min":round(f1s[1],3),"f1_max":round(f1s[2],3)},
 "ancoramento_cdr3_media":{"id_denovo_puro":round(mr[0],3),"id_ancorado_banco":round(ma[0],3),
     "ganho":round(ma[0]-mr[0],3)},
 "ms_peptideos_media":round(nms[0],1),
 "por_cadeia_semente0":detailed["rows"],
}
json.dump(report, open(os.path.join(OUT,"poc_resultado.json"),"w"), indent=2, ensure_ascii=False)

print(f"Banco: {len(DB)} cadeias ({sum(1 for v in DB.values() if v['source']=='RNAseq_cand')} candidatos RNA-seq, "
      f"2 F59, {sum(1 for v in DB.values() if v['source']=='distrator')} distratores)")
print(f"MS simulada: ~{nms[0]:.0f} peptídeos de novo/rodada · {len(ENZ)} proteases · "
      f"cobertura {P_OBS:.0%} · erro de novo {P_SUB:.0%} · {N_SEEDS} sementes\n")
print("=== (1) MS confirma candidatos expressos (semente 0) ===")
print(f"{'cadeia':15}{'fonte':12}{'verdade':8}{'pept.disc':10}{'enz':4}{'CDR3conf':9}{'CHAMADA':8}")
for r in sorted(detailed["rows"], key=lambda r:(-r['call'],r['id'])):
    print(f"{r['id']:15}{r['source']:12}{'EXPR' if r['truth'] else '-':8}{r['discrim']:^10}{r['enz']:^4}"
          f"{('sim('+str(r['cdr3_nenz'])+'e)') if r['cdr3_conf'] else '-':9}{'EXPRESSO' if r['call'] else '-':8}")
print(f"\nMédia em {N_SEEDS} sementes:  Precisão={precs[0]:.2f}  Recall={recs_[0]:.2f}  "
      f"F1={f1s[0]:.2f} (min {f1s[1]:.2f} – max {f1s[2]:.2f})")
print("\n=== (2) Banco corrige CDR3 ruidosa da MS ===")
print(f"identidade média à CDR3 verdadeira:  de novo puro = {mr[0]:.2f}  ->  ancorado no banco = {ma[0]:.2f}"
      f"   (ganho +{ma[0]-mr[0]:.2f})")
print(f"\nSaídas -> {OUT}/poc_resultado.json  ·  {OUT}/ms_sim/")
