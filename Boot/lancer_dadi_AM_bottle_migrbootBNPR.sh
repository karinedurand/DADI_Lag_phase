#!/bin/bash
#SBATCH -p workq

module load bioinfo/dadi/1.6.3_modif

    python ./script_inference_demo_new_models_folded.2024BNPR.py -o AM_bottle_migr_BNPR${SLURM_ARRAY_TASK_ID} -y pop1 -x pop2  -p 120,130,140 -f thinned.BN_PR.vcf.gz.recode.vcf.data.boot${SLURM_ARRAY_TASK_ID}.txt.boot2024.fs -m AM_bottle_migr -l -v 


