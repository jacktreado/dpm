#!/bin/bash

#SBATCH -p scavenge
#SBATCH -t 0-01:00:00
#SBATCH --job-name=a2j_process
#SBATCH -o out/a2j_process.out
#SBATCH --mem-per-cpu=10G

# load MATLAB
module load MATLAB

# user's netid
netid="$USER"

# matlab input information
NCELLS=64
nsmall=24
calA0=$1
kb=0
trun=$2
T0=$3

# sim info string
jamstr=a2j_N"$NCELLS"_n"$nsmall"_calA0"$calA0"_kb"$kb"_trun"$trun"_T0"$T0"

# simulation directory string
simloc=/gpfs/loomis/pi/ohern/"$netid"/dpm/jam
simstr=$simloc/$jamstr/"$jamstr"_seed

# matlab string
savedir=$simloc/matfiles
mkdir -p $savedir
savestr=$savedir/"$jamstr".mat


# str to give to matlab
MCODE="addpath ~/dpm/viz/jam/; processJammedDPMEnsemble('$simstr','$savestr'); quit"

echo running in matlab
matlab -nodisplay -r "$MCODE"
echo finished running in matlab


