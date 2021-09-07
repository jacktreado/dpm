#!/bin/bash

# directories with code
gitdir=~/dpm
srcdir=$gitdir/src
maindir=$gitdir/main/meso2D

# directories for dpm + tumor simulations
dpmdir=/gpfs/loomis/project/fas/ohern/"$USER"/dpm
tumordir=$dpmdir/meso2D

# parent directory for simulation output
outputdir=$dpmdir/mesoPenta2D

# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out

# inputs about simulation
n1=$1
calA0=$2
dh=$3
kb0=$4
betaEff=$5
ctcdel=$6
ctch=$7
cL=$8
aL="${9}"
cB="${10}"
cKb="${11}"

# inputs about cluster
partition="${12}"
time="${13}"
startSeed="${14}"
endSeed="${15}"

# other params
dh=1e-3

# name strings
basestr=penta_n"$n1"_ca"$calA0"_kb0"$kb0"_be"$betaEff"_cd"$ctcdel"_ch"$ctch"_cL"$cL"_aL"$aL"_cB"$cB"_cKb"$cKb"
runstr="$basestr"_startseed"$startSeed"_endseed"$endSeed"

# make directory specific for this simulation
simdatadir=$outputdir/$basestr
mkdir -p $simdatadir

# create paramf
paramf="$simdatadir"/"$basestr".params
echo n1=$n1 >> $paramf
echo calA0=$calA0 >> $paramf
echo dh=$dh >> $paramf
echo kb0=$kb0 >> $paramf
echo betaEff=$betaEff >> $paramf
echo ctcdel=$ctcdel >> $paramf
echo ctch=$ctch >> $paramf
echo cL=$cL >> $paramf
echo aL=$aL >> $paramf
echo cB=$cB >> $paramf
echo cKb=$cKb >> $paramf
echo partition=$partition >> $paramf
echo time=$time >> $paramf
echo startSeed=$startSeed >> $paramf
echo endSeed=$endSeed >> $paramf
cat $paramf



# compile into binary using packing.h
binf=bin/"$runstr".o
mainf=$maindir/mesoPenta2D.cpp

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

# LOOP OVER SEEDS
for seed in `seq $startSeed $endSeed`; do
    # print to console
    echo Adding seed = $seed / $endSeed to task file $taskf

    # create output file
    positionFile=$simdatadir/"$basestr"_seed"$seed".pos

    # create runString
    runString="./$binf $n1 $calA0 $dh $kb0 $betaEff $ctcdel $ctch $cL $aL $cB $cKb $seed $positionFile"

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
echo \#SBATCH -t $time >> $slurmf
echo sed -n \"\$\{SLURM_ARRAY_TASK_ID\}p\" "$taskf" \| /bin/bash >> $slurmf
cat $slurmf

# run sbatch file
echo -- running on slurm in partition $partition
sbatch $slurmf


# ====================
#       INPUTS
# ====================
# # inputs about simulation
# n1
# calA0
# dh
# kb0
# betaEff
# ctcdel
# ctch
# cL
# aL
# cB
# cKb
# partition
# time
# startSeed
# endSeed










