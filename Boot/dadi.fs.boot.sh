#!/bin/bash

#SBATCH -p workq
#SBATCH --mem=50G


module load bioinfo/dadi/1.6.3_modif

# Répertoire contenant les fichiers

#perl /work/kdurand/dgimi/DADI/VCF/convert_vcf_to_dadi_input.pl thinned.FL_I.vcf.gz.recode.vcf /work/kdurand/dgimi/DADI/VCF/INFL.pop
python /work/user/kdurand/dgimi/DADI/boot2024/INPR_boot.py


#perl /work/kdurand/dgimi/DADI/VCF/convert_vcf_to_dadi_input.pl thinned.BR_I.vcf.gz.recode.vcf /work/kdurand/dgimi/DADI/VCF/INBR.pop
#python /work/kdurand/dgimi/DADI/VCF/INBR.py

#perl /work/kdurand/dgimi/DADI/VCF/convert_vcf_to_dadi_input.pl thinned.MS_I.vcf.gz.recode.vcf  /work/kdurand/dgimi/DADI/VCF/INMS.pop
#python /work/kdurand/dgimi/DADI/VCF/INMS.py


#perl /work/kdurand/dgimi/DADI/VCF/convert_vcf_to_dadi_input.pl thinned.PR_I.vcf.gz.recode.vcf /work/kdurand/dgimi/DADI/VCF/INPR.pop
#python /work/kdurand/dgimi/DADI/VCF/INPR.py



