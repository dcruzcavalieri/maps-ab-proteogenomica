# (a) Dataset público pareado RNA-seq/BCR-seq + Espectrometria de Massas — VALIDADO

**Objetivo:** dataset público com, do **mesmo doador**, (i) sequenciamento de repertório
de células B (BCR-seq/RNA-seq → candidatos de sequência) **e** (ii) MS de IgG de
soro/plasma — para a validação proteogenômica (candidatos viram banco de busca da MS;
a MS confirma quais são expressos).

> **Status:** ✅ **validado online** (NCBI SRA + EBI PRIDE + Europe PMC), 2 out 2026.

---

## ★ ESCOLHIDO — dataset pareado completo (mesmos doadores)

**Park J, Bartzoka F, … Lavinder JJ, … Ippolito GC, Georgiou G.**
"Molecular features of the serological IgG repertoire elicited by egg-based,
cell-based, or recombinant haemagglutinin-based seasonal influenza vaccines."
***Lancet Microbe*, 2024.** PMID 39667375 · PMC11807745 · DOI 10.1016/j.lanmic.2024.06.002.
(grupo Georgiou/Ippolito, UT Austin — pioneiros do "Ig-seq")

| Camada | Acessão (confirmada) | Conteúdo |
|---|---|---|
| **Repertório BCR (transcriptômica)** | **NCBI BioProject `PRJNA1023018`** — **105 runs SRA** | FASTQ brutos de **Illumina MiSeq + PacBio** (PacBio = anticorpo **full-length**). |
| **MS de IgG sérico** | **PRIDE `PXD046473`** | Bottom-up LC-MS/MS de IgG purificado por afinidade a HA (eluato/flow-through, dia 0 / dia 28). |
| **Verdade-terreno (mAbs)** | GenBank (VH/VL) · PDB `8UZC` · EMD `EMD-42839` | Monoclonais validados e estruturas de alguns clones. |

- **15 doadores** em 3 coortes: A1–5 (RIV4, HA recombinante), B1–5 (eIIV4, ovo), C1–5 (ccIIV4, célula).
- **Por que é ideal:** mesmo doador com **BCR-seq (inclusive PacBio full-length → candidatos de domínio V completo)** + **MS de IgG sérico** + **monoclonais de sequência conhecida** como verdade-terreno. É exatamente o desenho proteogenômico que propusemos — e o PacBio dá candidatos full-length (melhor que o short-read da amostra COVID original).

---

## ★ COMPANHEIRO — benchmark de MS com verdade-terreno (lado proteômico)

**Chernigovskaya M, et al.** "Benchmarking of mass spectrometry-based antibody
repertoire profiling reveals extensive systematic biases." ***Cell Systems*, 2025**
(DOI 10.1016/j.cels.2025.101449 — **já citado no TCC do Gabriel**).

- **PRIDE `PXD055846`** — **6 anticorpos monoclonais de sequência conhecida**, em 70 combinações de concentração/policlonalidade (com e sem fundo de IgG policlonal), **4 proteases** (tripsina, quimotripsina, AspN, quimo+tripsina), **3.360 corridas LC-MS/MS**; análise com MaxQuant/MSFragger + **de novo (Casanovo)**.
- **Uso:** verdade-terreno para o **benchmark de acurácia da WS1** e para **calibrar o lado MS** e o modelo de erro de novo da POC in-silico (b). Alinha com o MAPS-ab (multi-protease, de novo, fundo policlonal).

---

## Como acessar (para o pipeline / equipe)
- **SRA (PRJNA1023018):** NCBI SRA Run Selector → `prefetch` + `fasterq-dump` (ou ENA FASTQ). 105 runs; MiSeq (short) + PacBio (full-length).
- **PRIDE (PXD046473 / PXD055846):** página do projeto em https://www.ebi.ac.uk/pride/archive/projects/PXD046473 (e …/PXD055846) → FTP dos `.raw`/`.mzML`.

## Observação de escopo
O dataset pareado é de **influenza HA** (não arbovírus/COVID). Ele é o **padrão-ouro
para desenvolver e validar o método** proteogenômico. Depois de validado, o mesmo
método se aplica à amostra COVID original e, no futuro, aos dados de arbovírus do
guarda-chuva do MAPS-ab.

## Verificações realizadas nesta etapa
- Europe PMC: paper, open access (PMC11807745), databanks PDB+PRIDE, data-sharing com PRJNA1023018.
- NCBI SRA: `PRJNA1023018` → **105 runs** (confirma BCR-seq MiSeq+PacBio disponível).
- PRIDE: `PXD046473` (pareado, B cell + plasma) e `PXD055846` (benchmark, 6 mAbs conhecidos) confirmados.
