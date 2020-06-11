if [ -e msdss_new_adjust_v2.sh ]
then
rm msdss_new_adjust_v2.sh
fi

touch msdss_new_adjust_v2.sh
echo "jid_01=\$(swarm -f /data/barbourcr/somalogic_msdss_modeling/scripts/msdss_new_adjust_v2/swarmR_iter0.sw -t 16 -g 12 --partition=quick,norm --job-name msdss_new_adjust_v2 --module R)" >> msdss_new_adjust_v2.sh
echo "jid_02=\$(sbatch --depend=afterok:\$jid_01 --mem=6g --partition=quick,norm --job-name=msdss_new_adjust_v2 /data/barbourcr/somalogic_msdss_modeling/scripts/msdss_new_adjust_v2/results_batch_iter0)" >> msdss_new_adjust_v2.sh

for j in {1..119}
do

i=$((j-1))

echo "jid_${j}1=\$(swarm -f /data/barbourcr/somalogic_msdss_modeling/scripts/msdss_new_adjust_v2/swarmR_iter${j}.sw -t 16 -g 12 --partition=quick,norm --depend=afterok:\$jid_${i}2 --job-name msdss_new_adjust_v2 --module R)" >> msdss_new_adjust_v2.sh
echo "jid_${j}2=\$(sbatch --depend=afterok:\$jid_${j}1 --mem=6g --partition=quick,norm --job-name=msdss_new_adjust_v2 /data/barbourcr/somalogic_msdss_modeling/scripts/msdss_new_adjust_v2/results_batch_iter${j})" >> msdss_new_adjust_v2.sh
done