//////////////////////////////
// Strong convergence test  //
//////////////////////////////
//
// First, the script generates the reference Brownian process with time step dt.
// Then, it computes N trajectories using the stochastic splitting method 
//  at N different time steps DT >= dt.
// For each time step, the stochastic forces are found 
//  from the reference Brownian process.
// The trajectories are saved to files.
//
// The stochastic pendulum and PQ splitting are used in this example script.

#include "methods_particle.h"

int main(int argc, char* argv[])
{	
	lfloat q0 = strtoflt128(argv[1], NULL);
	lfloat p0 = strtoflt128(argv[2], NULL);
	TEMP = strtoflt128(argv[3], NULL);
	lfloat eta = strtoflt128(argv[4], NULL);
	lfloat time = strtoflt128(argv[5], NULL);
	int split_order = std::stoi(argv[7]);
	int stoch_order = std::stoi(argv[8]);
	int N_SAMPLES = std::stoi(argv[9]);
	int seed = std::stoi(argv[10]);
	int N = std::stoi(argv[11]);
	int Num_of_steps[N];
	for(int i = 0; i < N; i++) {
		Num_of_steps[i] = std::stoi(argv[12+i]);
	}	
	bool strong_conv = true;
	lfloat dt = time/Num_of_steps[0];

///////////////////////////////////////////////////// 	

	Evolve sol(q0, p0, eta, dt, argv[6], strong_conv, split_order, stoch_order);
	sol.generate_Brownian_process(Num_of_steps[0], dt, seed);

	for( int i=0; i<N; i++ ) {
		std::vector<lfloat> qt;
		std::vector<lfloat> pt;
		DT = time/Num_of_steps[i];
		sol.vars[0] = q0;
		sol.vars[1] = p0;
		sol.evolve(Num_of_steps[i], N_SAMPLES, qt, pt);
		save_data(qt, "Out/Pendulum/SelfConv/PQ/pendulum_qt", Num_of_steps[i], split_order, stoch_order, seed);
		save_data(pt, "Out/Pendulum/SelfConv/PQ/pendulum_pt", Num_of_steps[i], split_order, stoch_order, seed);
	}
	return 0;
}
