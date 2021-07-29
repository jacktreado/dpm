#!/bin/bash

# directories with code
gitdir=~/dpm
srcdir=$gitdir/src
maindir=$gitdir/main/tumor2D

# directories for dpm + tumor simulations
dpmdir=/gpfs/loomis/project/fas/ohern/"$USER"/dpm
tumordir=$dpmdir/tumor2D

# parent directory for simulation input, throw error if it does not exist
inputdir=$tumordir/intInit
if [[ ! -d $inputdir ]]; then
    echo ERROR: directory $inputdir does not exist, so ending.
    exit 1
fi

# parent directory for simulation output
outputdir=$dpmdir/intInvade

# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out

# inputs about initialization file
aN=$1
aCalA0=$2
tCalA0=$3
areaRatio=$4

# inputs about invasion
l1=$5
l2=$6
v0=$7
Dr=$8
Ds=$9
kecm="${10}"
ecmbreak="${11}"
dDr="${12}"
dPsi="${13}"
Drmin="${14}"

# inputs about cluster
partition="${15}"
time="${16}"
startSeed="${17}"
endSeed="${18}"

# name strings
inputstr=intInit_aN"$aN"_ac"$aCalA0"_tc"$tCalA0"_aR"$areaRatio"
basestr=intInvade_aN"$aN"_ac"$aCalA0"_tc"$tCalA0"_aR"$areaRatio"_l1"$l1"_l2"$l2"_v0"$v0"_Dr"$Dr"_dPsi"$dPsi"
runstr="$basestr"_startseed"$startSeed"_endseed"$endSeed"

# make directory specific for this simulation
simdatadir=$outputdir/$basestr
mkdir -p $simdatadir

# compile into binary using packing.h
binf=bin/"$runstr".o
mainf=$maindir/interfaceInvasion.cpp

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

# get list of files
flist="$inputdir"/"$inputstr"_seed*.pos
let arrsz=0
for f in $flist; do
    let arrsz=$arrsz+1
done
if [[ $arrsz -lt 2 ]]
then
    echo flist = $flist
    echo arrsz = $arrsz, which is too small, ending.
    exit 1
else
    echo flist has $arrsz files, adding to task list...
fi

# create task file
taskf=tasks/"$runstr".task
rm -f $taskf

# loop over files
let fcount=0

# LOOP OVER FILES. 
for f in flist; do
    # parse file name
    file=${f##*/}
    baseid=${file%%.pos}
    seed=${baseid#*seed*}
    echo seed = $seed, input file = $file

    # check if seed is in correct range
    if [[ $seed -lt $startSeed ]]
    then
        echo seed = $seed too small, skipping...
        continue
    else
        # check if file is empty
        if [[ ! -s $f ]]
        then
            echo file $file is empty, skipping...
            continue
        else
            # increment file count
            let fcount=$fcount+1
            echo seed = $seed is ready for primetime, adding to task file.
        fi
    fi

    # create output file
    outputf=$simdatadir/"$basestr"_seed"$seed".vdos

    # create runString
    runString="./$binf $f $kl $kb $vdosf"

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










