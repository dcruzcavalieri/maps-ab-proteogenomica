#!/usr/bin/env bash
# ============================================================================
# Download do BCR-seq pareado — BioProject PRJNA1023018 (Park/Georgiou 2024)
# Executar na infraestrutura da equipe (precisa de sra-tools; dados ~dezenas de GB).
# NÃO roda no sandbox do Claude (rede restrita + volume grande).
# ============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
DEST="${1:-$HERE/dados/bcrseq_raw}"
mkdir -p "$DEST"

# ---- PRIORIDADE: 15 PacBio full-length (um por doador A1-5/B1-5/C1-5) ----
# São a melhor fonte de CANDIDATOS de domínio V completo (VH/VL full-length).
PACBIO="$HERE/dados/bcrseq_pacbio_runs.txt"
echo ">> PacBio full-length (15 runs) -> $DEST"
while read -r acc; do
  [ -z "$acc" ] && continue
  echo "  prefetch $acc"; prefetch "$acc" -O "$DEST"
  fasterq-dump "$DEST/$acc/$acc.sra" -O "$DEST" --split-files --skip-technical
done < "$PACBIO"

# ---- OPCIONAL: 90 Illumina MiSeq (amplicon VH/VL pareado, maior profundidade) ----
# Descomente para baixar também o short-read (TRUST4 de alta cobertura).
# ILL="$HERE/dados/bcrseq_illumina_runs.txt"
# while read -r acc; do prefetch "$acc" -O "$DEST"; \
#   fasterq-dump "$DEST/$acc/$acc.sra" -O "$DEST" --split-files --skip-technical; done < "$ILL"

echo ">> Pronto. FASTQ em $DEST"
echo "   Próximo: reconstruir candidatos (TRUST4 nos MiSeq; para PacBio, usar os"
echo "   reads full-length direto como candidatos VH/VL)."
