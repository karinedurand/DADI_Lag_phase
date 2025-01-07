#!/bin/bash
#SBATCH -p workq
module load bioinfo/dadi/1.6.3_modif
python ./script_inference_demo_new_models_folded.2024.py -o SI_${SLURM_ARRAY_TASK_ID} -y pop1 -x pop2 -p 120,130,140 -f INPR2024.fs -m SI -l -v  

#python ./script_inference_demo_new_models_folded.py -o SI_id -y pop1 -x pop2 -p 80,90,100  -f pop1_pop2_NP_78.fs -m SI -l -v  &>> SI_id.log
#python ./script_inference_demo_new_models_folded.py -o SC_id -y pop1 -x pop2 -p 80,90,100  -f pop1_pop2_NP_78.fs -m SI -l -v  &>> SC_id.log
