#!/bin/bash

# directories with code
gitdir=~/dpm
srcdir=$gitdir/src
maindir=$gitdir/main/meso2D

# directories for dpm + tumor simulations
dpmdir=/gpfs/loomis/project/fas/ohern/"$USER"/dpm
outputdir=$dpmdir/mesoPenta2D

# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p bin
mkdir -p slurm
mkdir -p out

# inputs about simulation
n1=$1
calA0=$2
kb0=$3
betaEff=$4
ctcdel=$5
ctch=$6
cL=$7
aL=$8
cB=$9
cKb="${10}"

# inputs about cluster
partition="${11}"
time="${12}"

# redo compilation
redocomp="${13}"

# other params
dh=1e-3
seed=1

# name strings
basestr=penta_n"$n1"_ca"$calA0"_kb0"$kb0"_be"$betaEff"_cd"$ctcdel"_ch"$ctch"_cL"$cL"_aL"$aL"_cB"$cB"_cKb"$cKb"
runstr="$basestr"_run

# compile into binary using packing.h
binf=bin/mesoPenta2D.o
mainf=$maindir/mesoPenta2D.cpp

# run compiler
if [[ $redocomp == 1 || ! -f $binf ]]
then
    echo compiling with : g++ --std=c++11 -O3 -I "$srcdir" "$mainf" "$srcdir"/*.cpp -o $binf
    g++ --std=c++11 -O3 -I "$srcdir" "$mainf" "$srcdir"/*.cpp -o $binf 
else
    echo $binf already exists, running on slurm
fi

# check compilation
if [[ ! -f $binf ]]
then
    echo -- binary file does not exist, compilation failed.
    exit 1
fi

# create output file
positionFile=$outputdir/"$basestr".pos
bondFile=$outputdir/"$basestr".bnd

# create runString
runString="./$binf $n1 $calA0 $dh $kb0 $betaEff $ctcdel $ctch $cL $aL $cB $cKb $seed $positionFile $bondFile"

# setup slurm files
slurmf=slurm/"$runstr".slurm
job_name="$runstr"
runout=out/"$runstr".out
rm -f $slurmf

# echo about time
echo -- running time = $time for $partition

echo -- PRINTING SLURM FILE...
echo \#\!/bin/bash >> $slurmf
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo \#SBATCH -o $runout >> $slurmf
echo \#SBATCH -t $time >> $slurmf
echo $runString >> $slurmf
cat $slurmf

# run sbatch file
echo -- running on slurm in partition $partition
sbatch $slurmf


# ====================
#       INPUTS
# ====================
# n1
# calA0
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










