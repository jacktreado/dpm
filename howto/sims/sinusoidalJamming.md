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