# The Deformable Particle Model
By Jack Treado, Yale University

Welcome to the Deformable Particle Model (DPM)!

This repository broadly allows for the simulation of deformable particles. Currently (as of date: "`r format(Sys.time(), '%d %B, %Y')`") we only support 2D particles, but 3D is on the way!



References that use the DPM:

J. D. Treado, D. Wang, A. Boromand, M. P. Murrell, M. D. Shattuck, and C. S. O'Hern, "Bridging particle deformability and collective response in soft solids," _Phys. Rev. Materials_ **5** 055605 (2021).

A. Boromand, A. Signoriello, J. Lowensohn, C. S. Orellana, E. R. Weeks, F. Ye, M. D. Shattuck, and C. S. O'Hern, "The role of deformability in determining the structural and mechanical properties of bubbles and emulsions," _Soft Matter_ **15** 5854 (2019).

A. Boromand, A. Signoriello, F. Ye, C. S. O'Hern, and M. D. Shattuck, "Jamming of deformable polygons," _Phys. Rev. Lett._ **121** 248003 (2018).


## Compilation

Once the repository has been downloaded, any code can be compiled using the `g++` compiler with the `--std=c++11` flag.

Specific simulations are stored in the `main` directory, which all use files written in the `src` directory.

To compile a given simulation to a binary (say `bin.o`), make the `dpm` directory your working directory and use:

`g++ -O3 --std=c++11 -I src main/[DIR NAME]/[MAIN FILE NAME].cpp src/*.cpp -o bin.o`



# How-to guides

## Simulations

### Jamming of deformable particles with sinusoidal preferred angles

Simulations with particles of arbitrarily-lobed, deformable particles. See [sinusoidalJamming.md](howto/sims/sinusoidalJamming.md)










