#!/bin/bash
#SBATCH --account=def-sibiryak
#SBATCH --array=0-999
#SBATCH --time=23:59:59
#SBATCH --mem=4000
#SBATCH --cpus-per-task=1
#SBATCH --job-name=fvd
#SBATCH --output=./Out/slurm.out

###########################
#  Decay rate measurement #
###########################
#
# The script runs a series of real-time stochastic simulations of the real massive scalar field 
#      with quartic self-interaction (unstable potential). 
# The field is prepared in the state of thermal equilibrium around the false vacuum.
# It is evolved until the decay is detected (or simulation times out).
# The decay times are recorded; the decay rate is determined by measuring the slope of the survival probability.
#
# This is a slurm script that was run on a cluster of Digital Research Alliance of Canada (alliancecan.ca).
# Command to compile the cpp-script:
# g++ -std=c++17 -o decay_rate.out decay_rate.cpp methods_fields.cpp -lfftw3 -lm

SIZE=100                       # lattice size
N_x=8192                       # number of lattice points
DT=0.1                         # time step
N_ENS=10                       # number of simulations per process
TIMESPAN=4000                  # maximum simulation time

BETA=4                         # splitting order
GAMMA=3                        # strong order

TEMP=0.1                       # temperature
ETA=10                         # damping coefficient

declare -a VARS=($SIZE $N_x $DT $N_ENS $TIMESPAN $BETA $GAMMA $TEMP $ETA)

./decay_rate.out "${VARS[@]}"
