#!/bin/bash
#SBATCH -p workq

output_file="parameters_INPR_SI.csv"
#SI    nu1, nu2, Ts = params
echo "X nu1 nu2 b g1 m12 m21 migr Tam Tb Ts Tsc likelihood" > "$output_file"

for ((X=1; X<=300; X++)); do
  folder="SI_${X}_2024_*"
  file_path="/work/user/kdurand/dgimi/DADI/2024_new_repet_3models/INPR/$folder/*.txt"
  
  # Vérifier si le dossier contient des fichiers correspondants
  files_exist=($file_path)
  if [[ -e "${files_exist[0]}" ]]; then
    echo "Processing folder: $folder"
    
    # Extraire le bloc de la section "BFGS"
    bfgs_block=$(grep -A 6 "Optimization : 'BFGS'" $file_path)
    bfgs_line=$(echo $bfgs_block | sed 's/"\n"/, /g')
    
    # Extraire les paramètres et la likelihood optimisés 
    optimized_params=$(echo "$bfgs_line" | awk -F'[][]' '{ print $2 }' | tr -d ',')
    optimized_likelihood=$(echo "$bfgs_block" | grep "Optimized log-likelihood" | awk '{print $NF}')
    
    # Extraire les paramètres individuels
    nu1=$(echo "$optimized_params" | awk '{print $1}')
    nu2=$(echo "$optimized_params" | awk '{print $2}')
    ts=$(echo "$optimized_params" | awk '{print $3}')
    
    # Écrire les paramètres et la likelihood dans le fichier CSV
    echo   "SI_INPR_$X $nu1 $nu2 $b $g1 $migr $tam $tb $ts $tsc $optimized_likelihood" >> "$output_file"
  fi
done

echo "Extraction des paramètres terminée. Le fichier $output_file a été créé."
