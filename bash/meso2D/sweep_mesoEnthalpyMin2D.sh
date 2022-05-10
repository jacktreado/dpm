#!/bin/bash

# directory for all output for cell simulations
outputdir=/gpfs/loomis/project/fas/ohern/jdt45/dpm

# directory for simulations specific to sphereGel
simtypedir=$outputdir/mesoEnthalpyMin2D

# directory to save matfiles
savedir=$simtypedir/matfiles

# inputs
forceRecalc=$1
partition=$2
time=$3

# get list of all data directories 
flist=($(ls -d "$simtypedir"/*))

# also get list of matfiles
mflist=($(ls "$savedir"/*_processed.mat))

# Count number of directories, error if empty
let nf=0
for f in "${flist[@]}"; do
	let nf=$nf+1
done
if [[ $nf -eq 0 ]]; then
	echo Error, number of files = $nf, ending here.
	exit 1
else
	if [[ $forceRecalc -eq 1 ]]; then
		echo Found $nf files in flist, now computing matfiles for all sim dirs.
	else
		echo Found $nf files in flist, now finding ones that dont have an associated matfile.
	fi
fi

# Count number of matfiles, set forceRecalc to 1 if empty
let nmf=0
for f in "${mflist[@]}"; do
	let nmf=$nmf+1
done
if [[ $nmf -eq 0 && $forceRecalc -eq 0 ]]; then
	echo No matfiles found, so setting forceRecalc to 1, now sweeping through all simulations
fi

# loop over directories, find matfile name, call process file
if [[ $forceRecalc -eq 1 ]]; then
	echo Looping over all $nmf files, will recalc all matfiles, including $nmf extant matfiles.
else
	echo Looping over all $nmf files, but will skip $nmf extant matfiles.
fi

# create task file
taskf=tasks/sweep_mesoEnthalpyMin2D.task
rm -f $taskf

let nfsubmit=0
for f in "${flist[@]}"; do
	# get sim name & details
	simdatadir=$f
	simname=${f##*/}
	searchstr="$simname"_seed

	# skip if you have found matfiles directory
	if [[ $simdatadir == $savedir ]]; then
		echo Found matfiles directory, skipping.
		continue
	fi

	# tmp matfile
	mftmp="$simname""_processed.mat"

	# if not recalc, check if matfile already exists
	if [[ $forceRecalc != 1 ]]; then
		mffound=0
		for mf in "${mflist[@]}"; do
			# get test matfile name
			mftest=${mf##*/}

			# test matfile
			if [[ $mftest == $mftmp ]]; then
				echo "\+\+" Found matfile $mftest, skipping
				mffound=1
				break
			fi
		done

		# continue loop if found
		if [[ $mffound == "1" ]]; then
			continue
		else
			echo Processing sims to $mftmp.
		fi
	else
		echo Processing sims to $mftmp.
	fi

	# add to task file
	let nfsubmit=$nfsubmit+1
	MCODE="addpath ~/dpm/viz/meso2D; processMesoEnthalpyMin2D('$simdatadir','$searchstr','$mftmp'); quit"
	echo matlab -nodisplay -r \""$MCODE"\" >> $taskf
done

# test if task file was created
if [[ ! -f "$taskf" ]]
then
    echo task file not created, ending before job submission
    exit 1
fi

# setup slurm files
slurmf=slurm/sweep_mesoEnthalpyMin2D.slurm
job_name=sweep_mesoEnthalpyMin2D
runout=out/sweep_mesoEnthalpyMin2D-%a.out
rm -f $slurmf

# echo about time
echo -- running time = $time for $partition

echo -- PRINTING SLURM FILE...
echo \#\!/bin/bash >> $slurmf
echo \#SBATCH --cpus-per-task=1 >> $slurmf
echo \#SBATCH --array=1-$nfsubmit >> $slurmf
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo \#SBATCH -o $runout >> $slurmf
echo module load MATLAB >> $slurmf
echo sed -n \"\$\{SLURM_ARRAY_TASK_ID\}p\" "$taskf" \| /bin/bash >> $slurmf
cat $slurmf

# run sbatch file
echo -- running on slurm in partition $partition
sbatch -t $time $slurmf














