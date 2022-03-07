#!/bin/bash

# directory for all output for cell simulations
outputdir=/gpfs/loomis/project/fas/ohern/jdt45/dpm

# directory for simulations specific to sphereGel
simtypedir=$outputdir/mesoBoxGrowth


# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p $simtypedir
mkdir -p tasks

# simulation inputs
N=$1
NGROWTH=$2
calA0_base=$3
kb_base=$4
dp0=$5
da0=$6
cB=$7
bbreak=$8
th0_min_scale=$9
P0="${10}"

# name strings
basestr=mesoBoxGrowth_N"$N"_NG"$NGROWTH"_ca"$calA0_base"_kb"$kb_base"_dp"$dp0"_da"$da0"_cB"$cB"_bb"$bbreak"_t0m"$th0_min_scale"_P0"$P0"

# get mafile string to save data
savestr="$simtypedir"/"$basestr".mat
taskf=tasks/mesoBoxGrowth.task

# create matlab command
MCODE="addpath ~/dpm/viz/meso2D; mesoBoxGrowth('$savestr',$N,$NGROWTH,$calA0_base,$kb_base,$dp0,$da0,$cB,$bbreak,$th0_min_scale,$P0); quit"

# add to task file
if [[ ! -f "$taskf" ]]
then
	echo Creating new task file $taskf, adding line $MCODE
	echo module load MATLAB/2021a";" matlab -nodisplay -r \""$MCODE"\" >> $taskf
else
	echo Task file $taskf exists, appending line $MCODE
	echo module load MATLAB/2021a";" matlab -nodisplay -r \""$MCODE"\" >> $taskf
fi



# ====================
#       INPUTS
# ====================
# 
# 1. N
# 2. NGROWTH
# 3. calA0
# 4. kb
# 5. dp0
# 6. da0
# 7. cB
# 8. bbreak
# 9. th0_min_scale (sets min th0 in units of 0.5 * pi)
# 10. P0

# NOTE: this only adds to task file mesoBoxGrowth.task in tasks/ directory
# ** If running a new batch file, delete task file (rm tasks/mesoBoxGrowth.task)
# ** Else, jobs will append original task file

# !!!! MAKE SURE TO CHECK IF YOU WANT TO DELETE TASK FILE !!!!!
#
#
#
#
#




