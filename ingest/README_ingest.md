# (ii) Ingest dos dados reais pareados — runbook

Objetivo: levar o dataset proteogenômico pareado validado na etapa (a) até a
**integração real** (run_integration.py), produzindo (1) quais candidatos de
anticorpo estão expressos e (2) as CDR3 ancoradas/corrigidas.

> O **download bruto e o TRUST4 rodam na infraestrutura da equipe** (volume grande,
> rede). Os scripts abaixo foram escritos para isso. A ferramenta de integração já
> foi **validada ponta a ponta** (na MS simulada da POC) e está pronta.

## Datasets (confirmados na etapa a)
| Camada | Acesso | Volume | Observação |
|---|---|---|---|
| BCR-seq (candidatos) | NCBI **PRJNA1023018** — 105 runs | dezenas de GB | **15 PacBio full-length** (1/doador A1–5,B1–5,C1–5) + 90 Illumina MiSeq (amplicon VH/VL). |
| MS de IgG sérico | PRIDE **PXD046473** | ~92 GB total | 29 `.raw` (~36 GB) **+ 32 `*_psms.txt`** (peptídeos já identificados — usar estes). FTP: `ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD046473`. |
| Benchmark (gabarito MS) | PRIDE **PXD055846** | ~839 GB | 6 mAbs de sequência conhecida, 4 proteases, de novo (Casanovo). Só os resultados de busca/de novo p/ calibração. |

## Arquivos gerados nesta etapa
- `dados/bcrseq_PRJNA1023018_manifest.tsv` — 105 runs (run, plataforma, modelo, estratégia, doador).
- `dados/bcrseq_pacbio_runs.txt` (15) · `dados/bcrseq_illumina_runs.txt` (90) — listas p/ download.
- `ingest/download_bcrseq.sh` · `ingest/download_ms.sh` — scripts de download.
- `ingest/airr_to_candidates.py` — AIRR (TRUST4/pipeline) → FASTA de candidatos.
- `ingest/prep_ms_peptides.py` — `_psms.txt`/MaxQuant/PEAKS → CSV (Peptide, ALC, Area).
- `ingest/run_integration.py` — **a integração** (candidatos + MS → expressos + CDR3 ancorada). ✅ validada.

## Receita ponta a ponta (por doador)
```bash
cd MAPS-ab-proteogenomica
# 1) BCR-seq (prioridade PacBio full-length)
bash ingest/download_bcrseq.sh dados/bcrseq_raw

# 2) Candidatos de sequência (VH/VL):
#    (a) PacBio: reads full-length já são candidatos → traduzir p/ FASTA; ou
#    (b) TRUST4 nos MiSeq → AIRR → FASTA:
python ingest/airr_to_candidates.py resultado_trust4_airr.tsv -o dados/candidatos_A1.fasta --min-consensus 2

# 3) MS (baixar as tabelas de PSM e normalizar, nomeando por protease)
bash ingest/download_ms.sh dados/ms_raw
python ingest/prep_ms_peptides.py dados/ms_raw/PXD046473/7649_JL_5a_psms.txt -o dados/ms_pep/A1/A1_S1_Elu_Trypsin.csv
#   (repita por amostra/protease; inclua a protease no nome do arquivo)

# 4) INTEGRAÇÃO (por doador): candidatos + MS -> expressos + CDR3 ancorada
python ingest/run_integration.py --candidates dados/candidatos_A1.fasta \
       --ms-dir dados/ms_pep/A1 --min-alc 70 --out docs/integ_A1.json
```

## Pareamento por doador
- PacBio (candidatos): runs rotulados por doador (A1–5 RIV4, B1–5 eIIV4, C1–5 ccIIV4).
- MS (PXD046473): amostras rotuladas por doador + dia (S1=dia0, S2=dia28) + fração
  (Elu = eluato anti-HA, FT = flow-through). Casar doador a doador.

## Validação
- **Gabarito interno:** os monoclonais validados do estudo (VH/VL em GenBank) devem
  aparecer como **expressos** e com CDR3 batendo — teste de acurácia direto.
- **Calibração do lado MS:** rodar o benchmark `PXD055846` (6 mAbs conhecidos) pela
  mesma `run_integration.py` para medir sensibilidade por protease/concentração.

## Status
- Metadados e manifests: ✅ obtidos e salvos.
- Scripts de download + conversores + integração: ✅ escritos; integração ✅ validada
  na MS simulada da POC (4/4 expressos corretos, 7/7 rejeitados).
- Pendente (infra da equipe): executar download + TRUST4 + normalização e rodar a
  integração nos dados reais.
