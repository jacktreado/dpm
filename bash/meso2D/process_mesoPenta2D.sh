#!/bin/bash

# inputs
partition=$1
time=$2

# slurm file name
basestr=process_mesoPenta2D
slurmf=slurm/"$basestr".slurm
job_name="$basestr"
runout=out/"$basestr".out

# input info
floc=/gpfs/loomis/project/fas/ohern/jdt/dpm/mesoPenta2D
fpattern=mesoPenta2D
savestr="$floc"/"$basestr".mat

# matlab code
MCODE="addpath ~/dpm/viz/meso2D; processMesoPenta2D('$floc','$fpattern','$savestr'); quit"

# echo about time
echo -- running time = $time for $partition

echo -- PRINTING SLURM FILE...
echo \#\!/bin/bash >> $slurmf
echo \#SBATCH --cpus-per-task=1 >> $slurmf
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo \#SBATCH -o $runout >> $slurmf
echo \#SBATCH -t $time >> $slurmf
echo \#SBATCH --mem-per-cpu=50G >> $slurmf
echo module load MATLAB >> $slurmf
echo matlab -nodisplay -r \""$MCODE"\" >> $slurmf
cat $slurmf

# run sbatch file
echo -- running on slurm in partition $partition
sbatch $slurmf