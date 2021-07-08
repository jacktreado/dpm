#!/bin/bash
# directories with code

#example call: bash bash/epi2D/submit_laserAblation.sh 24 24 1.08 0.7 0.9 1.0 0.5 0.5 1.0 1000 pi_ohern,day,scavenge 0-12:00:00 1 1
cellsdir=~/dpm
srcdir=$cellsdir/src
maindir=$cellsdir/main/epi2D

# directory for all output for cell simulations
outputdir=/gpfs/loomis/project/fas/ohern/at965/dpm

# directory for simulations specific to laserAblation
simtypedir=$outputdir/ablate

# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p $simtypedir
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out

# inputs
# NT = simulation time (in tau), time = human time
NCELLS=$1
NV=$2
calA0=$3
phiMin=$4
phiMax=$5
kl=$6
att=$7
B=$8
Dr0=$9
NT="${10}"
partition="${11}"
time="${12}"
numRuns="${13}"
startSeed="${14}"

numSeedsPerRun=1
let numSeeds=$numSeedsPerRun*$numRuns
let endSeed=$startSeed+$numSeeds-1

# name strings
basestr=ablate_N"$NCELLS"_NV"$NV"_calA0"$calA0"_kl"$kl"_att"$att"_B"$B"_Dr0$Dr0
runstr="$basestr"_startseed"$startSeed"_endseed"$endSeed"

# make directory specific for this simulation
simdatadir=$simtypedir/$basestr
mkdir -p $simdatadir

# compile into binary using packing.h
binf=bin/"$runstr".o
mainf=$maindir/laserAblation.cpp

echo Running laserAblation simulations with parameters:
echo NCELLS = "$NCELLS"
echo NV = "$NV"
echo calA0 = "$calA0"
echo phiMin = "$phiMin"
echo phiMax = "$phiMax"
echo kl = "$kl"
echo att = "$att"
echo B = "$B"
echo Dr0 = "$Dr0"

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
        energyf=$simdatadir/$filestr.energy

        # append to runString
        runString="$runString ; ./$binf $NCELLS $NV $calA0 $phiMin $phiMax $kl $att $B $Dr0 $runseed $NT $posf $energyf"
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
echo \#SBATCH --mail-type=END,FAIL >> $slurmf
echo \#SBATCH --mail-user=andrewtondata@gmail.com >> $slurmf
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
# 2. NV
# 3. calA0
# 4. phiMin
# 5. Ptol
# 6. kl
# 7. kb
# 8. att
# 9. number of timesteps
# 10. partition
# 11. time
# 12. number of runs (number of array entries, i.e. arraynum)
# 13. start seed (end seed determined by number of runs)
