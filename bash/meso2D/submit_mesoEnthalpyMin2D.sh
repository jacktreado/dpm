#!/bin/bash

# directories with code
cellsdir=~/dpm
srcdir=$cellsdir/src
maindir=$cellsdir/main/meso2D

# directory for all output for cell simulations
outputdir=/gpfs/loomis/project/fas/ohern/jdt45/dpm

# directory for simulations specific to jamming
simtypedir=$outputdir/mesoEnthalpyMin2D

# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p $simtypedir
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out

# inputs
NCELLS=$1
n1=$2
calA0=$3
kb0=$4
betaEff=$5
da0=$6
dl0=$7
P0=$8
ctch=$9
cL="${10}"
cB="${11}"
partition="${12}"
time="${13}"
numRuns="${14}"
startSeed="${15}"

# other variables
numSeedsPerRun=1

let numSeeds=$numSeedsPerRun*$numRuns
let endSeed=$startSeed+$numSeeds-1

# name strings
basestr=mesoHMin2D_N"$NCELLS"_n"$n1"_ca"$calA0"_kb0"$kb0"_be"$betaEff"_da"$da0"_dl"$dl0"_P"$P0"_h"$ctch"_cL"$cL"_cB"$cB"
runstr="$basestr"_startseed"$startSeed"_endseed"$endSeed"

# make directory specific for this simulation
simdatadir=$simtypedir/$basestr
mkdir -p $simdatadir

# compile into binary using packing.h
binf=bin/"$runstr".o
mainf=$maindir/mesoEnthalpyMin2D.cpp
echo Running mesoEnthalpyMin2D simulations with parameters:
echo NCELLS = "$NCELLS"
echo n1 = "$n1"
echo calA0 = "$calA0"
echo kb0 = "$kb0"
echo betaEff = "$betaEff"
echo da0 = "$da0"
echo dl0 = "$dl0"
echo P0 = "$P0"
echo ctch = "$ctch"
echo cL = "$cL"
echo cB = "$cB"

# run compiler
rm -f $binf
g++ --std=c++11 -O3 -I "$srcdir" "$mainf" "$srcdir"/*.cpp -o $binf 
echo compiling with : g++ --std=c++11 -O3 -I "$srcdir" "$mainf" "$srcdir"/*.cpp -o $binf   

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
for seed in `seq $startSeed $numSeedsPerRun $endSeed`; do
    # count files
    let fcount=$fcount+1

    # echo to console
    echo On base seed $seed

    # echo string of numSeedPerRun commands to task file
    runString="cd `pwd`"

    # loop over seeds to go into runString
    let ssMax=$numSeedsPerRun-1

    for ss in `seq 0 $ssMax`; do
        # get seed for actual run
        let runseed=$seed+ss

        # get file str
        filestr="$basestr"_seed"$seed"

        # create output files
        posf=$simdatadir/$filestr.posctc

        # append to runString
        runString="$runString ; ./$binf $NCELLS $n1 $calA0 $kb0 $betaEff $da0 $dl0 $P0 $ctch $cL $cB $runseed $posf"
    done

    # finish off run string
    runString="$runString ;"

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
# 1. NCELLS
# 2. n
# 3. calA0
# 4. kb0
# 5. betaEff
# 6. da0
# 7. dl0
# 8. P0
# 9. ctch
# 10. cL (perimeter aging)
# 11. cB (bending angle aging)
# 12. partition
# 13. time
# 14. number of runs (number of array entries, i.e. arraynum)
# 15. start seed (end seed determined by number of runs)










