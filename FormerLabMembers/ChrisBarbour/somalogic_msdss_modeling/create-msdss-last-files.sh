cd /data/barbourcr/somalogic_msdss_modeling/scripts/msdss_last_new_adjust_v2/

for j in {0..119}
do

if [ -e msdss_last_new_adjust_v2_ratios_forest_iter${j}_.R ]
then
rm msdss_last_new_adjust_v2_ratios_forest_iter${j}_.R
fi

cp msdss_last_new_adjust_v2_ratios_forest.R msdss_last_new_adjust_v2_ratios_forest_iter${j}_.R

echo y | Rswarm -rfile=msdss_last_new_adjust_v2_ratios_forest_iter${j}_.R --sfile=seeds.txt --path=/data/barbourcr/somalogic_msdss_modeling/scripts/msdss_last_new_adjust_v2/ --start=0 --reps=10 --sims=1 --ext1=.rds &> /dev/null
rm msdss_last_new_adjust_v2_ratios_forest_iter${j}_.sw

if [ -e swarmR_iter${j}.sw ]
then
rm swarmR_iter${j}.sw 
fi

touch swarmR_iter${j}.sw

for i in {1..10}
do
echo "R --vanilla --args var_path \"./scripts/msdss_last_new_adjust_v2/iter${j}_ratios.txt\" remove_outliers FALSE < /data/barbourcr/somalogic_msdss_modeling/scripts/msdss_last_new_adjust_v2/msdss_last_new_adjust_v2_ratios_forest_iter${j}_${i}.R" >> swarmR_iter${j}.sw
done

if [ -e results_batch_iter${j} ]
then
rm results_batch_iter${j} 
fi

i=$((j+1))

touch results_batch_iter${j}

echo "#!/bin/bash
# This file is runR
# 
#SBATCH -J runR
date
module load R
R --vanilla --args remove_outliers FALSE path1 "./scripts/msdss_last_new_adjust_v2/msdss_last_new_adjust_v2_ratios_forest_iter${j}_" pull_iter "./scripts/msdss_last_new_adjust_v2/iter${j}" push_iter "./scripts/msdss_last_new_adjust_v2/iter${i}" < /data/barbourcr/somalogic_msdss_modeling/scripts/msdss_last_new_adjust_v2/msdss_last_new_adjust_v2_ratios_forest_results.R" >> results_batch_iter${j}

done