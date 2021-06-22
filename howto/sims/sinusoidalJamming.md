# Jamming of deformable particles with sinusoidal preferred angles

Using the main file `main/jam/bidisperseSinusoidalParticleJamming.cpp`, you can generate jammed packings (at a specified pressure) of bidisperse, purely repulsive, deformable particles with lobed shapes sets by a sinusoidally-varying preferred angle profile. 

See some example jammed configurations below.

<p float="left">
  <img src="/img/jammedTriangles.png" width="400" />
  <img src="/img/jammedPentalobes.png" width="400" /> 
</p>

Particles drawn with empty circles are "rattlers" i.e. are not part of the force-bearing backbone of the packing. 

## Input

Once you compile the main file, run the simulation using the generated binary file. The input parameters are:
* `NCELLS`: integer number of particles
* `nsmall`: integer number of vertices on small particles, `n` on large particles is set by the 1.4:1.0 size ratio. 
* `calA0`: **preferred** shape parameter of all particles
	* defined as `calA0 = p_0^2/(4 * pi * a_0)`, where `p_0` and `a_0` are the preferred perimeter and areas of the particles, respectively
	* Note that input `calA0` should be greater than or equal to 1 always. 
* `kl`: mechanical constant for perimeter
* `kb`: mechanical constant for curvature
* `thA`: amplitude of sinusoidal preferred angles
* `thK`: wavenumber of sinusoidal preferred angles (**sets number of lobes**)
* `Ptol`: pressure tolerance, jammed configs with final pressure `P` satisfy `Ptol  < P < 2 * Ptol`
* `Ftol`: force tolerance of energy minimization protocol
* `seed`: integer to seed random number generator
* `positionFile`: path to file to store position data for jammed configuration. 

To generate the above images, run a compiled binary `bin.o` using

* _left_: `./bin.o 12 24 1.04 1.0 0.01 3.0 3.0 1e-7 1e-12 1 pos.test`

* _right_: `./bin.o 12 24 1.20 1.0 0.01 10.0 5.0 1e-7 1e-12 1 pos.test`

## Output

The final input to the main file, `positionFile` (above called `pos.test`), will store the positions of all vertices in the jammed state as well as some other useful information about the state of the system. 

Each row in this file starts with a five-character keyword that denotes what information that row contains. The keywords are:
* `NEWFR`: starts frame
	* If the code prints multiple configurations, this can be used to identify the start of a new configuration.
* `NUMCL`: number of cells
* `PACKF`: packing fraction (particle areas / box area)
* `BOXSZ`: box lengths `Lx` and `Ly`
* `STRSS`: components of the virial stress tensor, sorted as `Sxx`, `Syy`, `Sxy`
* `CINFO`: information about cell
* `VINFO`: information about vertex 
* `ENDFR`: ends frame

### CINFO 

There is one `CINFO` row per cell. For a given cell `mu`, the row is formatted as:

`CINFO` `nv` `zc` `zv` `a0` `a` `p`

where
* `nv` number of vertices on cell `mu`
* `zc` number of cell-cell contacts on cell `mu`. 
* `zv` number of vertex-vertex contacts on cell `mu`.
* `a0` preferred area of cell `mu`
* `a` instantaneous area of cell `mu`
* `p` instantaneous perimeter of cell `mu`

### VINFO

There is also one `VINFO` row per vertex. For a given vertex `i` on cell `mu`, the row is formatted as:

`CINFO` `mu` `i` `x` `y` `r` `l0` `t0`

where `mu` and `i` refer to the cell and vertex indices, respectively, and
* `x`: x coordinate of vertex
* `y`: y coorindate of vertex
* `r`: radius of vertex
* `l0`: preferred length of edge between vertices `i` and `i+1`
* `t0`: preferred angle for angle centered on vertex `i`





# Running on the cluster

To run an ensemble of jammed configurations on a remote cluster, use the following steps. **Note**: currently only supports clusters that use the [Slurm workload manager](https://slurm.schedmd.com/). 

### **READ THIS FIRST**

BEFORE YOU SUBMIT ANYTHING, make sure to edit line `4` to set the `netid` variable to your remote cluster username. This will place all data in a folder named `~/project/dpm` if using the Yale clusters. If you have a different file storage system, edit the variable `outputdir` on line `12` to be your desired location for data output. 


## Using the `submit_` bash script

To clone the repository to your remote cluster, `ssh` onto the login node and use 
```bash 
>> git clone https://github.com/jacktreado/dpm.git
``` 

The bash script [submit_bidisperseSinusoidalParticleJamming.sh](/bash/jam/submit_bidisperseSinusoidalParticleJamming.sh) will submit an ensemble of simulations to the cluster as an array job. 

1. Once you have logged in to the remote cluster, `cd` to `~/dpm/bash/jam`.
2. Use `>> cat submit_bidisperseSinusoidalParticleJamming.sh` to see the input list, which is:
```bash
# ====================
#       INPUTS
# ====================
# 1. NCELLS
# 2. n
# 3. calA0
# 4. kb
# 5. thA (lobey preferred angle amplitude)
# 6. thK (lobey preferred angle wavenumber)
# 7. partition
# 8. time
# 9. number of runs (number of array entries, i.e. arraynum)
# 10. start seed (end seed determined by number of runs)
```
3. Jobs are submitted with `>> bash submit_bidisperseSinusoidalParticleJamming.sh [INPUTS]` where `[INPUTS]` must match the list above. See input description below.

## Bash script inputs

1. `NCELLS`: integer number of particles
2. `n`: integer number of vertices on small particles (`nsmall` above)
3. `calA0`: **preferred** shape parameter of all particles
4. `kb`: mechanical constant for curvature
5. `thA`: amplitude of sinusoidal preferred angles
6. `thK`: wavenumber of sinusoidal preferred angles (**sets number of lobes**)
7. `partition`: Slurm partition to queue + run jobs
8. `time`: Slurm time limit
9. `numSeedsPerRun`: number of jobs in batch submission
10. `startSeed`: first seed in array, last seed is `startSeed + numSeedsPerRun`

## Data storage

By default, data will be stored in `~/project/dpm/jam` *if you have set* `netid` *to your netid*. 

Each file will begin with `lobes_`, with the parameters for the ensemble following. For example, an ensemble with `NCELLS=16`, `n=24`, `calA0=1.10`, `kb0.01`, `thA=3.0`, and `thK=3.0` will be stored in the directory
```bash
~/project/dpm/jam/lobes_N16_n24_calA01.10_kb0.01_thA3.0_thK3.0
```
Each file will begin with the string `lobes_N16_n24_calA01.10_kb0.01_thA3.0_thK3.0` and will end with `_seed[seedNumber].pos`, where `seedNumber` denotes the seed used in the random number generator. This differentiates different members of the particular ensemble. 

## Data processing (Slurm scheduler only)

The MATLAB function [processJammedDPMEnsemble](/viz/jam/processJammedDPMEnsemble.m) will aggregate statistical data stored the ensemble folder (i.e. `~/project/dpm/jam/lobes_`) and save it to a `.mat` file (see [MATLAB documentation](https://www.mathworks.com/help/matlab/ref/matlab.io.matfile.html) stored in `~/project/dpm/jam/matfiles`. 

To submit this function to the cluster, use the Slurm script [slurm_bidisperseSinusoidalParticleJamming.slurm](/bash/jam/slurm_bidisperseSinusoidalParticleJamming.slurm). This script must be edited via the command line, so open the file using a command line editor like `vim` or `emacs`. 

**NOTE**: You MUST edit the `netid` variable to be your netid, or else the data won't be saved in the right location!

Slurm-specific options, like `partition` or `time` are controlled by the `#SBATCH` headers. See [this Slurm cheatsheet](https://slurm.schedmd.com/pdfs/summary.pdf) for all cluster options. 

Simulation options are set by these variables:
```bash
# matlab input information
NCELLS=16
nsmall=24
calA0=1.10
kb=1.0
thA=3.0
thK=3.0
```

To process a given set of simulations, edit these variables & the Slurm headers appropriately, then use `sbatch` to submit the script to the cluster:
```bash
>> sbatch slurm_bidisperseSinusoidalParticleJamming.slurm
```

A `matfile` with ensemble statistics will be saved in the folder `~/project/dpm/jam/matfiles` with the same name as the ensemble folder.