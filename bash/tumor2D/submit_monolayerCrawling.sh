#!/bin/bash

# directories with code
gitdir=~/dpm
srcdir=$gitdir/src
maindir=$gitdir/main/tumor2D

# directories for dpm + tumor simulations
dpmdir=/gpfs/loomis/project/fas/ohern/"$USER"/dpm
tumordir=$dpmdir/tumor2D
outputdir=$tumordir/mono

# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out

# inputs
NT=$1
tN=$2
tNV=$3
tCalA0=$4
l1=$5
Dr=$6
gamtt=$7
tau=$8
partition=$9
time="${10}"
numSeeds="${11}"
startSeed="${12}"

# other variables
NPRINTSKIP=5e3

# determine total number of seeds
let endSeed=$startSeed+$numSeeds-1

# name strings
basestr=mono_NT"$NT"_N"$tN"_n"$tNV"_tc"$tCalA0"_l1"$l1"_Dr"$Dr"_gamtt"$gamtt"_tau"$tau"
runstr="$basestr"_startseed"$startSeed"_endseed"$endSeed"

# make directory specific for this simulation
simdatadir=$outputdir/$basestr
mkdir -p $simdatadir

# compile into binary
binf=bin/"$runstr".o
mainf=$maindir/monolayerCrawling.cpp
rm -f $binf
echo compiling with : g++ --std=c++11 -O3 -I "$srcdir" "$mainf" "$srcdir"/*.cpp -o $binf
g++ --std=c++11 -O3 -I "$srcdir" "$mainf" "$srcdir"/*.cpp -o $binf    

# check compilation
if [[ ! -f $binf ]]
then
    echo -- binary file does not exist, compilation failed.
    exit 1
fi

# create task file
taskf=tasks/"$runstr".task
rm -f $taskf

# loop over files
let fcount=0

# LOOP OVER FILES. 
for seed in `seq $startSeed 1 $endSeed`; do
    # count files
    let fcount=$fcount+1

    # echo to console
    echo On base seed $seed

    # get file str
    filestr="$basestr"_seed"$seed"

    # create output files
    posf=$simdatadir/$filestr.pos

    # create command for simulation run
    runString="cd `pwd` ; ./$binf $NT $NPRINTSKIP $tN $tNV $tCalA0 $l1 $Dr $gamtt $tau $seed $posf"

    # echo to task file
    echo "$runString" >> $taskf
done

# test if task file was created
if [[ ! -f "$taskf" ]]
then
    echo task file not created, ending before job submission
    exit 1
fi

# get number of jobs to submit to each array
let arraynum=$fcount
echo -- total number of array runs = $arraynum

# setup slurm files
slurmf=slurm/"$runstr".slurm
job_name="$runstr"
runout=out/"$runstr"-%a.out
rm -f $slurmf

# echo about time
echo -- running time = $time for $partition

echo -- PRINTING SLURM FILE...
echo \#\!/bin/bash >> $slurmf
echo \#SBATCH --cpus-per-task=1 >> $slurmf
echo \#SBATCH --array=1-$arraynum >> $slurmf
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo \#SBATCH -o $runout >> $slurmf
echo sed -n \"\$\{SLURM_ARRAY_TASK_ID\}p\" "$taskf" \| /bin/bash >> $slurmf
cat $slurmf

# run sbatch file
echo -- running on slurm in partition $partition
sbatch -t $time $slurmf




# ====================
#       INPUTS
# ====================
# 1. NT
# 2. tN
# 3. tNV
# 4. tCalA0
# 5. l1
# 6. Dr
# 7. gamtt
# 8. tau
# 9. partition
# 10. time
# 11. number of seeds
# 12. start seed










