#!/bin/bash

# directory for all output for cell simulations
outputdir=/home/jdt45/scratch60

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
kl_base=$4
kb_base=$5
kc=$6
dl0=$7
da0=$8
cB=$9
P0="${10}"

# name strings
basestr=mesoBoxGrowth_N"$N"_NG"$NGROWTH"_ca"$calA0_base"_kl"$kl_base"_kb"$kb_base"_kc"$kc"_dl"$dl0"_da"$da0"_cB"$cB"_P0"$P0"

# get mafile string to save data
savestr="$simtypedir"/"$basestr".mat
taskf=tasks/mesoBoxGrowth.task

# create matlab command
MCODE="addpath ~/dpm/viz/meso2D; mesoBoxGrowth('$savestr',$N,$NGROWTH,$calA0_base,$kl_base,$kb_base,$kc,$dl0,$da0,$cB,$P0); quit"

# add to task file
if [[ ! -f "$taskf" ]]
then
	echo Creating new task file $taskf, adding line $MCODE
	echo matlab -nodisplay -r \""$MCODE"\" >> $taskf
else
	echo Task file $taskf exists, appending line $MCODE
	echo matlab -nodisplay -r \""$MCODE"\" >> $taskf
fi



# ====================
#       INPUTS
# ====================
# 
# 1. N
# 2. NGROWTH
# 3. calA0
# 4. kl
# 5. kb
# 6. kc
# 7. dl0
# 8. da0
# 9. cB
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




