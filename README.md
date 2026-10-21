# MAPS-ab-proteogenômica

Integração **proteogenômica** (RNA-seq + Espectrometria de Massas) para
**reconstrução e validação de anticorpos** — o braço transcriptômico do projeto
MAPS-ab. Combina candidatos de sequência vindos de RNA-seq/BCR-seq (o que o
paciente *pode* produzir) com leituras *de novo* de MS (o que *está expresso* como
proteína circulante): cada modalidade valida a outra.

- A MS fornece o **ground truth de expressão** (quais candidatos são reais).
- O transcriptoma fornece a **sequência completa** que ancora e corrige a MS.

Reusa o motor do MAPS-ab (`maps_ab.cdr3_repertoire`): recuperação de CDR3 ancorada
no FR4 + corroboração cross-enzima.

---

## Estrutura

```
MAPS-ab-proteogenomica/
├── busca_datasets/        (a) dataset público pareado RNA-seq + MS (validado)
│   └── a_shortlist_datasets_pareados.md
├── poc_insilico/          (b) prova de conceito in-silico da integração
│   ├── poc_integracao.py          # roda standalone (python poc_integracao.py)
│   ├── candidatos_pipeline.json    # 5 candidatos (dados embutidos)
│   ├── POC_relatorio.md
│   ├── poc_resultado.json          # métricas (geradas)
│   └── ms_sim/                     # peptídeos de novo simulados (gerados)
├── ingest/                (ii) ingest dos dados reais + integração
│   ├── run_integration.py          # a INTEGRAÇÃO (candidatos + MS -> resultado)
│   ├── airr_to_candidates.py       # AIRR (TRUST4) -> FASTA de candidatos
│   ├── prep_ms_peptides.py         # PSM/MaxQuant/PEAKS -> CSV (Peptide, ALC, Area)
│   ├── download_bcrseq.sh          # download SRA (infra da equipe)
│   ├── download_ms.sh              # download PRIDE (infra da equipe)
│   └── README_ingest.md            # runbook ponta a ponta
├── dados/                 manifests e listas de runs (dados brutos NÃO versionados)
├── docs/                  relatórios (.docx)
├── requirements.txt
└── README.md
```

---

## Requisitos

- **Python 3.10+** e **pandas** (`pip install -r requirements.txt`).
- **Pacote `maps_ab`** (fork *MAPS-ab-adaptado*) — **dependência externa**, não
  incluída aqui. Aponte para ele de uma das formas:
  - variável de ambiente: `export MAPS_AB=/caminho/para/MAPS-ab-adaptado`, ou
  - coloque o fork como **pasta irmã** deste repositório (`../MAPS-ab-adaptado`).
- Para a validação úmida (infra da equipe): `sra-tools` (prefetch/fasterq-dump),
  `TRUST4`, e acesso de rede ao NCBI/EBI.

---

## Início rápido — prova de conceito (etapa b)

Roda sozinha (dados embutidos), demonstra a integração com resposta conhecida:

```bash
pip install -r requirements.txt
export MAPS_AB=/caminho/para/MAPS-ab-adaptado   # ou pasta irmã
python poc_insilico/poc_integracao.py
```

Saída esperada: classificação da expressão **Precisão/Recall/F1 = 1,00** (12 sementes)
e correção da CDR3 **de novo 0,8 → ancorado 1,0**. Detalhes em
`poc_insilico/POC_relatorio.md`.

---

## Uso em dados reais (etapa ii)

Dataset pareado validado: **PRJNA1023018** (BCR-seq) + **PXD046473** (MS), mesmos
doadores (Park et al., *Lancet Microbe* 2024). Receita completa em
`ingest/README_ingest.md`. Resumo:

```bash
# 1) BCR-seq (prioridade PacBio full-length) — infra da equipe
bash ingest/download_bcrseq.sh dados/bcrseq_raw

# 2) Candidatos de sequência (VH/VL): TRUST4 -> AIRR -> FASTA
python ingest/airr_to_candidates.py resultado_trust4_airr.tsv -o dados/candidatos_A1.fasta --min-consensus 2

# 3) MS: baixar tabelas de PSM e normalizar (incluir a protease no nome do arquivo)
bash ingest/download_ms.sh dados/ms_raw
python ingest/prep_ms_peptides.py dados/ms_raw/PXD046473/7649_JL_5a_psms.txt -o dados/ms_pep/A1/A1_S1_Elu_Trypsin.csv

# 4) INTEGRAÇÃO por doador: candidatos + MS -> expressos + CDR3 ancorada
python ingest/run_integration.py --candidates dados/candidatos_A1.fasta \
       --ms-dir dados/ms_pep/A1 --min-alc 70 --out docs/integ_A1.json
```

### `run_integration.py` — entradas e saídas
- **Entrada:** `--candidates` (FASTA VH/VL) e `--ms-dir`/`--ms-csv` (peptídeos por
  protease; a protease é inferida do nome do arquivo — inclua `Trypsin`,
  `Chymotrypsin`, `LysC`, `AspN` ou `GluC`). Coluna de sequência aceita: `Peptide`/`Sequence`.
- **Saída (JSON):** por candidato — nº de peptídeos discriminantes, enzimas,
  CDR3 confirmada e a chamada **expresso/não**; além das CDR3 da MS **ancoradas** ao banco.

---

## Fontes de dados

| Camada | Acesso | Observação |
|---|---|---|
| BCR-seq (candidatos) | NCBI `PRJNA1023018` | 15 PacBio full-length + 90 Illumina MiSeq |
| MS de IgG sérico | PRIDE `PXD046473` | 29 `.raw` + 32 `_psms.txt` |
| Benchmark MS (gabarito) | PRIDE `PXD055846` | 6 mAbs conhecidos, 4 proteases, de novo |

> Os dados brutos (FASTQ/`.raw`, dezenas–centenas de GB) **não** são versionados;
> baixe-os com os scripts em `ingest/`.

---

## Limitações
- A POC é **in-silico** (MS simulada a partir de sequências conhecidas). A validação
  real depende de rodar a integração nos dados pareados acima (infra da equipe).
- O pareamento HC/LC por co-variação cross-paciente (MAPS-ab) requer coorte e é
  direção futura.

## Referências
- Park J. et al. *Lancet Microbe* 2024 (Ig-seq influenza HA) — PRJNA1023018, PXD046473.
- Chernigovskaya M. et al. *Cell Systems* 2025 (benchmark MS) — PXD055846.
- Bondt A. et al. *Cell Systems* 2021 (plasma IgG de novo; origem do F59).
- Song L. et al. *Nature Methods* 2021 — TRUST4.

## Licença
Sugerida: MIT (ver `LICENSE`). Ajuste conforme a política do seu grupo/UFES antes de publicar.
