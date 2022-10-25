#!/usr/bin/env sh

# python esm/scripts/extract.py esm2_t33_650M_UR50D test.fasta \
#   test_protein_esm2 --repr_layers 0 32 33 --include mean per_tok

# python esm/scripts/extract.py esm2_t33_650M_UR50D train.fasta \
#   train_protein_esm2  --include mean per_tok

python esm/scripts/extract.py esm2_t33_650M_UR50D thermo.fasta \
  thermo_protein_esm2  --include mean per_tok
