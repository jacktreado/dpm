#!/bin/bash

# directory for all output for cell simulations
outputdir=/gpfs/loomis/project/fas/ohern/jdt45/dpm

# directory for simulations specific to sphereGel
simtypedir=$outputdir/mesoDimerGrowth


# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p $simtypedir
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out

# simulation inputs
NGROWTH=$1
kl=$2
kb=$3
kc=$4
dl=$5
P0=$6

# inputs about cluster
partition=$7
time=$8

# name strings
basestr=mesoDimerGrowth_NG"$NGROWTH"_kl"$kl"_kb"$kb"_kc"$kc"_dl"$del_l0"_P0"$P0"

# get mafile string to save data
savestr="$simtypedir"/"$basestr".mat

# create matlab command
MCODE="addpath ~/dpm/viz/meso2D; mesoDimerGrowth('$savestr',$NGROWTH,$kl,$kb,$kc,$del_l0,$P0); quit"

# setup slurm files
slurmf=slurm/"$basestr".slurm
job_name="$basestr"
runout=out/"$basestr".out
rm -f $slurmf

# echo about time
echo -- running time = $time for $partition

echo -- PRINTING SLURM FILE...
echo \#\!/bin/bash >> $slurmf
echo \#SBATCH --cpus-per-task=1 >> $slurmf
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo \#SBATCH -o $runout >> $slurmf
echo \#SBATCH --mem-per-cpu=10G >> $slurmf
echo module load MATLAB/2021a >> $slurmf
echo matlab -nodisplay -r \""$MCODE"\" >> $slurmf
cat $slurmf

# run sbatch file
echo -- running on slurm in partition $partition
sbatch -t $time $slurmf



# ====================
#       INPUTS
# ====================
# 
# 1. NGROWTH
# 2. kl
# 3. kb
# 4. kc
# 5. dl
# 6. P0
# 7. partition
# 8. time


