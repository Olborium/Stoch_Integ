/*
    The script runs N_ENS real-time stochastic simulations of the scalar field 
with negative quartic self-interaction (unstable potential).
    The field is evolved until the decay is detected (i.e. the absolute value of the field
exceeds the threshold=10) or the time limit is reached.
    The decay times are saved to a file.
*/

#include "methods_fields.h"

int main(int argc, char* argv[])
{
    SIZE = std::stod(argv[1]);
    N_x = std::stoi(argv[2]);
    DT = std::stod(argv[3]);
    int N_ens = std::stoi(argv[4]);
    double timespan = std::stod(argv[5]);
    int split_order = std::stoi(argv[6]);
    int stoch_order = std::stoi(argv[7]);
    TEMP = std::stod(argv[8]);
    double eta = std::stod(argv[9]);
    DX = SIZE/N_x;

///////////////////////////////////////////////////// 
    std::vector<double> decay_times;

    Evolve sol(eta, split_order, stoch_order, -1);  // initialize the solver

    for(int i = 0; i < N_ens; i++) {
        sol.prepare_state();                  // set up the thermal initial state
        sol.evolve(decay_times, timespan);    // evolve the state for t=timespan or until the decay
        }                                     // if no decay, the decay time is set to zero
    std::random_device rd;
    int seed = std::abs(rd());
    save_data(2, decay_times, "Out/times/times", DT, TEMP, eta, seed);    // save the decay times 

    return 0;
}
