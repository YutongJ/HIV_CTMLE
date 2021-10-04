#!/bin/bash

##################### Change these constants ##############################
path_read="/home/yjin85/TMLE_RF/Codes/20210222_New_Sims_3AA"  # load directory
path_save="/home/yjin85/TMLE_RF/Results/20210222_New_Sims_3AA"  # write directory

# performing 1000 simulations
#splitting the whole simulation to subgroup
group_num=10 # the number of each sub-group
groups=100 # the number of groups

# if path directory doesn't exist, make it
[ ! -d ${path_read}/out ] && mkdir ${path_read}/out
[ ! -d ${path_read}/err ] && mkdir ${path_read}/err

for rhos in 0.75
do

for sig in 5
do

for s in 500 1000
do

for g in 0.01
do

for splt in {1..${groups}}
do

echo "#!/bin/bash" >> script.out
echo "#SBATCH --job-name=S_${splt}" >> script.out
echo "#SBATCH --partition=benkeser" >> script.out
echo "#SBATCH --output=${path_read}/out/S_${splt}.out" >> script.out
echo "#SBATCH --error=${path_read}/err/S_${splt}.err" >> script.out
echo "#SBATCH --mail-user=yjin85@emory.edu." >> script.out
echo "module purge" >> script.out
echo "module load R/4.0.2" >> script.out
echo "R CMD BATCH --no-restore '--args path_read=${path_read} path_save=${path_save} rhos=${rhos} sig=${sig} s=${s} g=${g} splt=${splt} group_num=${group_num}' ${path_read}/Sim_Multilevel_RF_4lel_TMLE.R" >> script.out

chmod +x script.out
sbatch ./script.out
rm -rf script.out

done
done
done
done
done