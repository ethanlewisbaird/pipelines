#!/usr/bin/env bash
# rds2h5ad.sh — one-shot Seurat .rds -> AnnData .h5ad converter.
#
# Usage: ./rds2h5ad.sh <input.rds> <output.h5ad> [keep_dump_dir]
#
# Steps:
#   1. dump_seurat.R streams matrices + metadata to a temp dir (R side).
#   2. reconstruct_h5ad.py assembles the .h5ad via memory-mapped reads.
#
# The intermediate dump dir is removed on success unless a 3rd arg is given.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <input.rds> <output.h5ad> [keep_dump_dir]" >&2
  exit 1
fi

IN_RDS="$1"
OUT_H5AD="$2"
KEEP="${3:-}"

if [[ ! -f "$IN_RDS" ]]; then
  echo "Input not found: $IN_RDS" >&2
  exit 1
fi

DUMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/rds2h5ad.XXXXXX")"
cleanup() {
  if [[ -z "$KEEP" ]]; then
    rm -rf "$DUMP_DIR"
  else
    echo "Dump dir kept at: $DUMP_DIR" >&2
  fi
}
trap cleanup EXIT

echo "==> [1/2] Dumping matrices from $IN_RDS (R) ..." >&2
Rscript "$HERE/dump_seurat.R" "$IN_RDS" "$DUMP_DIR"

echo "==> [2/2] Reconstructing $OUT_H5AD (Python) ..." >&2
python3 "$HERE/reconstruct_h5ad.py" "$DUMP_DIR" "$OUT_H5AD"

echo "==> Success: $OUT_H5AD" >&2
