# Deformable Particle Model in C++
By Jack Treado

This repository broadly allows for the simulation of deformable particles using an OOP, C++ simulation framework. Currently we only support 2D particles, but 3D is on the way!


## Compilation

Once the repository has been downloaded, any code can be compiled using the `g++` compiler with the `--std=c++11` flag.

Specific simulations are stored in the `main` directory, which all use files written in the `src` directory.

To compile a given simulation to a binary (say `bin.o`), make the `dpm` directory your working directory and use:

`g++ -O3 --std=c++11 -I src main/[DIR NAME]/[MAIN FILE NAME].cpp src/*.cpp -o bin.o`



# How-to guides

## Simulations

### Jamming of deformable particles with sinusoidal preferred angles

Simulations with particles of arbitrarily-lobed, deformable particles. See [sinusoidalJamming.md](howto/sims/sinusoidalJamming.md)










