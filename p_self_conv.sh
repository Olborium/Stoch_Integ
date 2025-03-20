#!/bin/bash

###########################
# Strong convergence test #
###########################
#
# The script tests the rate of convergence of the high-order stochastic schemes.
# It generates M groups of solutions computed with different time steps DT.
# Within each group, N solutions are computed with different time steps DT 
#    using the same Brownian process computed with the smallest of DT.
# All solutions are computed using the same strong and splitting orders.
# The average distance between the solutions computed at different DT is measured as a function of DT.
#
# Command to compile the cpp-script:
# g++ -std=c++17 -fext-numeric-literals -o p_self_conv.out p_self_conv.cpp methods_particle.cpp -lquadmath -lm
# Note that the cpp-script uses quadmath library that guarantees quadruple-precision arithmetics.

Q0=2            # initial coordinate
P0=2            # initial momentum
TEMP=1          # temperature
ETA=10          # damping coefficient
SIMTIME=100     # time of simulation
SAMPLES=10      # number of sample points (q,p) of a trajectory which are saved to a file
                # the samples are uniformly distributed between 0 and SIMTIME

INTEG="PQ"      # type of splitting method: PQ or LN
BETA=4          # order of splitting
GAMMA=3         # strong order

M=10            # size of ensemble
N=17            # number of different time steps DT

# total number of steps, DT[i] = SIMTIME/Ns[i], i = 0,..,N-1
declare -a Ns=("5000000" "2500000" "1000000" "500000" "250000" "100000" "50000" "25000" "10000" "5000" "2500" "1000" "500" "250" "200" "150" "100")

printf "\n beta=%d, gamma=%d, SimNo " $BETA $GAMMA
for ((j=0; j<M; j++)); do
    printf "%d, " $j
    declare -a VARS=($Q0 $P0 $TEMP $ETA $SIMTIME $INTEG $BETA $GAMMA $SAMPLES $j $N "${Ns[@]}")
    ./p_self_conv.out "${VARS[@]}"
done
printf " Calculation complete, processing... "
declare -a VARS=($TEMP $ETA $SIMTIME $BETA $GAMMA $SAMPLES $INTEG $M $N "${Ns[@]}")
python p_self_conv.py "${VARS[@]}"