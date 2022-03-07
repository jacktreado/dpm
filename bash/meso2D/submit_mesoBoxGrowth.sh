#!/bin/bash

mkdir -p out
mkdir -p slurm

partition=$1
time=$2

# figure out how many tasks are in task file
taskf=tasks/mesoBoxGrowth.task
if [[ ! -f "$taskf" ]]
then
	echo In submit_mesoBoxGrowth.sh, but no task file found! Ending here.
	exit 1
fi

# use wc to count lines
line_count=$(wc -l < $taskf)
let arraynum="$line_count"

# setup slurm files
slurmf=slurm/mesoBoxGrowth.slurm
job_name=mesoBoxGrowth
runout=out/mesoBoxGrowth-%a.out
rm -f $slurmf

echo -- PRINTING SLURM FILE...
echo \#\!/bin/bash >> $slurmf
echo \#SBATCH --cpus-per-task=1 >> $slurmf
echo \#SBATCH --array=1-$arraynum >> $slurmf
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo \#SBATCH -o $runout >> $slurmf
echo module load MATLAB/2021a >> $slurmf
echo sed -n \"\$\{SLURM_ARRAY_TASK_ID\}p\" "$taskf" \| /bin/bash >> $slurmf
cat $slurmf

# run sbatch file
echo -- running on slurm in partition $partition
sbatch -t $time $slurmf




