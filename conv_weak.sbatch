#!/bin/bash
#SBATCH --account=def-sibiryak
#SBATCH --time=0-23:59
#SBATCH --mem=4000
#SBATCH --cpus-per-task=1000
#SBATCH --job-name=ps
#SBATCH --output=./Out/slurm.out

###########################
#  Weak convergence test  #
###########################
#
# The script runs a series of real-time stochastic simulations of the real massive scalar field 
#      with quartic self-interaction (stable potential). 
# The simulations are run at N different time steps using the stochastic scheme with strong order GAMMA and splitting order BETA.
# The thermodynamic properties of the field (temperature, power spectrum) are measured at N_SAMPLE points during the run
#      and compared with the theory predictions.
#
# This is a slurm script that was run on a cluster of Digital Research Alliance of Canada (alliancecan.ca).
# Command to compile the cpp-script:
# g++ -std=c++17 -o conv_weak.out conv_weak.cpp methods_fields.cpp -lfftw3 -lm

SIZE=100                       # lattice size
N_x=8192                       # number of lattice points

N_SAMPLE=500                   # number of samples

BETA=4                         # splitting order 
GAMMA=3 				       # strong order

TEMP=0.1                       # temperature
ETA=1                          # damping coefficient

N=14                           # number of different time steps
# time steps:
declare -a DTs=("0.005" "0.007" "0.01" "0.02" "0.03" "0.04" "0.05" "0.07" "0.1" "0.2" "0.3" "0.4" "0.5" "0.7")
# time to evolve between sample points:
declare -a Ts=(  "0.5"   "0.49" "0.5"   "0.5" "0.48" "0.48" "0.5"  "0.49"  "0.5" "0.6" "0.6" "0.8" "1.0" "1.4")

# generate the power spectra in 1000 independent runs (run 10 copies of the script to reach the total number of simulations 10000):
for ((i=0; i<1000; i++)); do
	for ((j=0; j<N; j++)); do ./conv_weak.out $SIZE $N_x ${DTs[$j]} $N_SAMPLE ${Ts[$j]} $BETA $GAMMA $TEMP $ETA; done &
done 
wait

###########################
#    Data processing      #
###########################
#
# measure the effective temperature of long and short Fourier modes of the field:
#
python conv_weak_temp.py $SIZE $N_x $N_SAMPLE $BETA $GAMMA $TEMP $ETA $N "${DTs[@]}" "${Ts[@]}"
#
# measure the field power spectrum:
#
python conv_weak_power_spectrum.py $SIZE $N_x $N_SAMPLE $BETA $GAMMA $TEMP $ETA $N "${DTs[@]}" "${Ts[@]}"