#!/usr/bin/env bash
# ============================================================================
# Download da MS pareada — PRIDE PXD046473 (IgG sérico, Park/Georgiou 2024)
# e, opcional, PXD055846 (benchmark Chernigovskaya 2025, 6 mAbs conhecidos).
# Executar na infraestrutura da equipe. NÃO roda no sandbox do Claude.
# ============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
DEST="${1:-$HERE/dados/ms_raw}"
mkdir -p "$DEST/PXD046473" "$DEST/PXD055846"

FTP046="https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD046473"

# ---- RECOMENDADO primeiro: as 32 tabelas de PSM (_psms.txt) ----
# Peptídeos JÁ identificados por amostra (~16-29 MB cada). Alimentam a integração
# diretamente, sem reprocessar os .raw (que somam ~36 GB).
echo ">> PXD046473: tabelas de PSM (_psms.txt)"
# lista os arquivos *_psms.txt do projeto e baixa só eles
wget -q -O - "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD046473/files?pageSize=300" \
 | grep -oE '"fileName":"[^"]*_psms\.txt"' | sed -E 's/.*:"(.*)"/\1/' | sort -u \
 | while read -r f; do echo "  $f"; wget -nc -P "$DEST/PXD046473" "$FTP046/$f"; done

# ---- OPCIONAL: os 29 arquivos .raw (reprocessar de novo com PEAKS/Casanovo) ----
# wget -nc -r -nd -np -A '*.raw' -P "$DEST/PXD046473" "$FTP046/"

# ---- OPCIONAL: benchmark PXD055846 (grande, ~839 GB no total) ----
# Use os resultados de busca/Casanovo para calibração; evite baixar todos os .raw.
# FTP058="https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/12/PXD055846"
# (ajuste ano/mês conforme a data de publicação do projeto)

echo ">> Pronto. PSM tables em $DEST/PXD046473"
echo "   Converta cada _psms.txt para o formato de entrada da integração"
echo "   (colunas Peptide, ALC(%)/score, Area) — ver ingest/prep_ms_peptides.py."
