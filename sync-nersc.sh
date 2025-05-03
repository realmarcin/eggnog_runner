#!/usr/bin/env bash
set -eox pipefail

rsync -avz \
  --exclude __pycache__ --exclude .snakemake --exclude .idea --exclude .git --exclude .mypy_cache --exclude venv \
  --exclude KE_phenotype.egg-info \
  -e "ssh -A -i ~/.ssh/nersc" \
  ../KG-Microbe cjneely@perlmutter-p1.nersc.gov:/global/cfs/cdirs/kbase/ke_prototype/cjneely/KG-MICROBE/
