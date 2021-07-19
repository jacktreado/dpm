#!/bin/bash

# directories with code
cellsdir=~/dpm
srcdir=$cellsdir/src
maindir=$cellsdir/main/tumor2D

# directory for all output for cell simulations
outputdir=/gpfs/loomis/project/fas/ohern/"$USER"/dpm

# directory for simulations specific to jamming
simtypedir=$outputdir/tumor2D

# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p $simtypedir
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out

# inputs
aN=$1
aNV=$2
tNV=$3
aDisp=$4
tDisp=$5
aCalA0=$6
tCalA0=$7
areaRatio=$8
prt=$9
partition="${10}"
time="${11}"
numRuns="${12}"
startSeed="${13}"

# other variables
numSeedsPerRun=1
let numSeeds=$numSeedsPerRun*$numRuns
let endSeed=$startSeed+$numSeeds-1

# name strings
basestr=intInit_N"$NCELLS"_n"$n1"_ca"$calA0"_be"$betaEff"_cL"$cL"_cB"$cB"_cKb"$cKb"
runstr="$basestr"_startseed"$startSeed"_endseed"$endSeed"

# make directory specific for this simulation
simdatadir=$simtypedir/$basestr
mkdir -p $simdatadir

# compile into binary using packing.h
binf=bin/"$runstr".o
mainf=$maindir/interfaceCompression.cpp

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
        posf=$simdatadir/$filestr.pos

        # append to runString
        runString="$runString ; ./$binf $aN $aNV $tNV $aDisp $tDisp $aCalA0 $tCalA0 $areaRatio $prt $runseed $posf"
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
# 1. aN
# 2. aNV
# 4. tNV
# 5. aDisp
# 5. tDisp
# 6. aCalA0
# 7. tCalA0
# 8. areaRatio
# 9. prt
# 10. partition
# 11. time
# 12. number of runs (number of array entries, i.e. arraynum)
# 13. start seed (end seed determined by number of runs)










