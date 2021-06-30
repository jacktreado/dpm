#!/bin/bash

#SBATCH --job-name=fireTest
#SBATCH --output=output/firetest.txt
#SBATCH --ntasks=1
#SBATCH -p pi_ohern,day,scavenge
#SBATCH --time=10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=andrewtondata@gmail.com

./firetest2.o
