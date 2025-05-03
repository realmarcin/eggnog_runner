#!/bin/sh
#SBATCH --nodes=5
#SBATCH --qos=regular
#SBATCH --time=12:00:00
#SBATCH --account=m4689
#SBATCH -C cpu

python /global/cfs/cdirs/kbase/ke_prototype/cjneely/KG-MICROBE/KG-Microbe/kg_microbe/utils/parsl_emapper_runner.py
