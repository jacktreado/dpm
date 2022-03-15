#!/bin/bash

# directory for all output for cell simulations
outputdir=/gpfs/loomis/project/fas/ohern/jdt45/dpm

# directory for simulations specific to sphereGel
simtypedir=$outputdir/mesoEnthalpyMin2D

# directory to save matfiles
savedir=$simtypedir/matfiles

# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p $simtypedir
mkdir -p $savedir
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out

# inputs about initialization file
NCELLS=$1
n1=$2
calA0=$3

# inputs about network formation simulation
kb0=$4
betaEff=$5
ctch=$6
da0=$7
dl0=$8
cL=$9
cB="${10}"
t0min="${11}"
P0="${12}"

# inputs about cluster
partition="${13}"
time="${14}"

# name strings
basestr=mesoHMin2D_N"$NCELLS"_n"$n1"_ca"$calA0"_kb"$kb0"_be"$betaEff"_h"$ctch"_da"$da0"_dl"$dl0"_cL"$cL"_cB"$cB"_t0m"$t0min"_P"$P0"
runstr="$basestr"_PROCESS
searchstr="$basestr"_seed

# access directory specific for this simulation
simdatadir=$simtypedir/$basestr
if [[ ! -d $simdatadir ]]
then
    echo -- sim directory "$simdatadir" does not exist, ending.
    exit 1
fi

# get mafile string to save data
savestr="$savedir"/"$basestr"_processed.mat

# create matlab command
MCODE="addpath ~/dpm/viz/meso2D; processMesoEnthalpyMin2D('$simdatadir','$searchstr','$savestr'); quit"

# setup slurm files
slurmf=slurm/"$runstr".slurm
job_name="$runstr"
runout=out/"$runstr".out
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
# 1. NCELLS
# 2. n
# 3. calA0
# 4. kb0
# 5. betaEff
# 6. ctch
# 7. da0
# 8. dl0
# 9. cL
# 10. cB
# 11. t0min
# 12. P0
# 13. partition
# 14. time


