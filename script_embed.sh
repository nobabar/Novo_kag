#!/usr/bin/env sh

# python esm/scripts/extract.py esm2_t33_650M_UR50D test.fasta \
#   test_protein_esm2  --include mean

python esm/scripts/extract.py esm2_t33_650M_UR50D train.fasta \
  train_protein_esm2  --include mean

python esm/scripts/extract.py esm2_t33_650M_UR50D wildtype.fasta \
  wildtype_protein_esm2  --include mean
