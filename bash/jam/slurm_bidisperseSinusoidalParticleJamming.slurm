#!/bin/bash

#SBATCH -p scavenge
#SBATCH -t 0-01:00:00
#SBATCH --job-name=sinejamming
#SBATCH -o out/sinejamming.out
#SBATCH --mem-per-cpu=2000

# load MATLAB
module load MATLAB

# user's netid
netid=jdt45

# matlab input information
NCELLS=16
nsmall=24
calA0=1.10
kb=0.01
thA=3.0
thK=3.0

# sim info string
jamstr=lobes_N"$NCELLS"_n"$nsmall"_calA0"$calA0"_kb"$kb"_thA"$thA"_thK"$thK"

# simulation directory string
simloc=/gpfs/project/fas/ohern/"$netid"/dpm/jam
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


