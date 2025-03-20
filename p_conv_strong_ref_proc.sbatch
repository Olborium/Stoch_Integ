#!/bin/bash
#SBATCH --partition=defq
#SBATCH --job-name=BrownG
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=23:59:59
#SBATCH --output=./Out/slurm.out

###########################
# Strong convergence test #
###########################
#
# The script generates M=10 reference Brownian processes with time step dt = SIMTIME / N
# and M associated reference trajectories computed using the 1st strong order 2nd splitting order scheme.
# These trajectories are then used to measure the accuracy of the high-order stochastic schemes.
#
# This is a slurm script for an HPC.
# Command to compile the cpp-script:
# g++ -std=c++17 -fext-numeric-literals -o p_conv_strong_ref_proc.out p_conv_strong_ref_proc.cpp methods_particle.cpp -lquadmath -lm
# Note that the cpp-script uses quadmath library that guarantees quadruple-precision arithmetics.

Q0=2				# initial coordinate
P0=2				# initial momentum
TEMP=1  			# temperature
ETA=10  			# damping coefficient
SIMTIME=100         # simulation time
SAMPLES=100         # number of sample points (q,p) of a trajectory which are saved to a file
                    # the samples are uniformly distributed between 0 and SIMTIME
N=100000000         # total number of steps, dt = SIMTIME / N
INTEG="PQ"          # type of splitting method: PQ or LN
BETA=2              # order of splitting
GAMMA=1             # strong order

./p_conv_strong_ref_proc.out $Q0 $P0 $TEMP $ETA $SIMTIME $INTEG $BETA $GAMMA $SAMPLES "0" $N >> ./Out/out.txt &
./p_conv_strong_ref_proc.out $Q0 $P0 $TEMP $ETA $SIMTIME $INTEG $BETA $GAMMA $SAMPLES "1" $N >> ./Out/out.txt &
./p_conv_strong_ref_proc.out $Q0 $P0 $TEMP $ETA $SIMTIME $INTEG $BETA $GAMMA $SAMPLES "2" $N >> ./Out/out.txt &
./p_conv_strong_ref_proc.out $Q0 $P0 $TEMP $ETA $SIMTIME $INTEG $BETA $GAMMA $SAMPLES "3" $N >> ./Out/out.txt &
./p_conv_strong_ref_proc.out $Q0 $P0 $TEMP $ETA $SIMTIME $INTEG $BETA $GAMMA $SAMPLES "4" $N >> ./Out/out.txt &
./p_conv_strong_ref_proc.out $Q0 $P0 $TEMP $ETA $SIMTIME $INTEG $BETA $GAMMA $SAMPLES "5" $N >> ./Out/out.txt &
./p_conv_strong_ref_proc.out $Q0 $P0 $TEMP $ETA $SIMTIME $INTEG $BETA $GAMMA $SAMPLES "6" $N >> ./Out/out.txt &
./p_conv_strong_ref_proc.out $Q0 $P0 $TEMP $ETA $SIMTIME $INTEG $BETA $GAMMA $SAMPLES "7" $N >> ./Out/out.txt &
./p_conv_strong_ref_proc.out $Q0 $P0 $TEMP $ETA $SIMTIME $INTEG $BETA $GAMMA $SAMPLES "8" $N >> ./Out/out.txt &
./p_conv_strong_ref_proc.out $Q0 $P0 $TEMP $ETA $SIMTIME $INTEG $BETA $GAMMA $SAMPLES "9" $N >> ./Out/out.txt &
wait