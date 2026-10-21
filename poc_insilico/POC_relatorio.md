# (b) Prova de conceito in-silico — integração proteogenômica (RNA-seq + MS)

**Objetivo.** Demonstrar, com **resposta conhecida**, a ideia central do braço
MAPS-ab × TCC: combinar candidatos de anticorpo vindos de **RNA-seq** com leituras
**de novo de MS** para (1) confirmar quais candidatos estão de fato **expressos** como
proteína e (2) usar o banco de RNA-seq para **corrigir** as CDR3 ruidosas da MS.

> Reuso direto do MAPS-ab: as leituras de MS simuladas são processadas pelo motor
> `maps_ab.cdr3_repertoire` (recuperação ancorada no FR4 + corroboração cross-enzima).

---

## Montagem (com verdade-terreno)

- **Banco (candidatos de RNA-seq):** 11 cadeias — os **5 candidatos reconstruídos pelo
  pipeline do Gabriel** (assemble181, assemble95, assemble193, assemble726, assemble4017)
  + **F59** (HC/LC) + **4 distratores** (anticorpos humanos reais, não relacionados: 2D22, EDE1 C8).
- **Expressas no soro (verdade-terreno):** apenas 4 — o par dominante do pipeline
  (assemble181_0 HC + assemble95_0 LC) e o F59 (HC/LC). As demais 7 cadeias **não**
  estão no soro (candidatos "só-transcrito" + distratores).
- **MS simulada:** digestão in-silico das cadeias expressas com **3 proteases**
  (tripsina, quimotripsina, AspN), **cobertura parcial (60%)**, **erro de de novo de 8%**,
  ambiguidade I/L e clivagem perdida → ~120 peptídeos "de novo"/rodada, no formato que o
  MAPS-ab consome. Repetido em **12 sementes**.

---

## Resultado (1) — MS confirma os candidatos expressos

| cadeia | fonte | verdade | pept. discriminantes | enzimas | CDR3 confirmada | chamada |
|---|---|---|---:|---:|---|---|
| F59_HC | F59 | EXPR | 24 | 3 | sim | **EXPRESSO** |
| F59_LC | F59 | EXPR | 9 | 3 | sim | **EXPRESSO** |
| assemble181_0 | RNA-seq | EXPR | 20 | 3 | — | **EXPRESSO** |
| assemble95_0 | RNA-seq | EXPR | 8 | 3 | — | **EXPRESSO** |
| 2D22_HC / 2D22_LC | distrator | — | 0 | 0 | — | — |
| EDE1C8_HC / EDE1C8_LC | distrator | — | 0 | 0 | — | — |
| assemble193_0 / 726_0 / 4017_0 | RNA-seq | — | 0 | 0 | — | — |

**Precisão = 1,00 · Recall = 1,00 · F1 = 1,00** (idêntico em todas as 12 sementes, min = max = 1,00).

A MS separa perfeitamente os candidatos **expressos** dos **só-transcrito** e dos
distratores, usando peptídeos **discriminantes** (que batem em um único candidato do
banco) e a corroboração da CDR3 por múltiplas proteases.

## Resultado (2) — o banco corrige a CDR3 ruidosa da MS

Identidade média das CDR3 recuperadas em relação à CDR3 verdadeira:

| | de novo puro | ancorado no banco RNA-seq |
|---|---:|---:|
| identidade CDR3 | **0,82** | **1,00** |

**Ganho de +0,18.** Com erro de de novo realista (8%), a leitura crua acerta ~82% da
CDR3; ancorá-la ao candidato de RNA-seq mais próximo recupera a sequência exata. É a
segunda direção da validação mútua: o transcriptoma **conserta** o proteoma.

---

## Leitura

Os dois resultados demonstram, de ponta a ponta e com resposta conhecida, que a
integração proteogenômica resolve as duas fragilidades do estudo original do TCC:
a MS fornece o **ground truth de expressão** que faltava (quais candidatos são reais),
e o banco de RNA-seq fornece a **sequência completa e limpa** que ancora a MS. O motor
do MAPS-ab é reusado sem modificação.

## Limitações e próximo passo

- É **in-silico**: a MS é simulada a partir de sequências conhecidas (modelo de erro,
  cobertura e proteases realistas, mas não espectros reais). O F59, por estar no banco
  como sua sequência exata, é um controle "fácil"; o teste informativo são os candidatos
  de RNA-seq (assemble181/95), corretamente confirmados, e os só-transcrito/distratores,
  corretamente rejeitados.
- **Próximo (validação úmida — WS3):** rodar a mesma integração sobre o **dataset real
  pareado** validado na etapa (a): MS `PXD046473` + BCR-seq `PRJNA1023018` (Park/Georgiou,
  *Lancet Microbe* 2024), e calibrar o lado MS com o benchmark `PXD055846` (Chernigovskaya
  2025, 6 mAbs de sequência conhecida).

## Reprodutibilidade

- Código: `poc_integracao.py` (determinístico, 12 sementes).
- Saídas: `poc_resultado.json` (métricas) e `ms_sim/MSsim_{Trypsin,Chymotrypsin,AspN}.csv`
  (peptídeos de novo simulados, semente 0).
- Executar: `python poc_integracao.py` (requer o pacote `maps_ab` do fork MAPS-ab-adaptado no path).
