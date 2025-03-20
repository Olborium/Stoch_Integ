//////////////////////////////
// Strong convergence test  //
//////////////////////////////
//
// The script generates the reference Brownian process with time step dt
// Then it computes a particle trajectory driven by the Brownian process with time step dt
//
// The Brownian process, its moments and the trajectory are saved to files
//
// The stochastic pendulum is used in this example script.

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
	bool strong_conv = true;
	lfloat dt = time/N;

///////////////////////////////////////////////////// 	

	Evolve sol(q0, p0, eta, dt, argv[6], strong_conv, split_order, stoch_order);
	sol.generate_Brownian_process(N, dt, seed);

	std::vector<lfloat> qt;
	std::vector<lfloat> pt;
	DT = time/N;
	sol.vars[0] = q0;
	sol.vars[1] = p0;
	sol.evolve(N, N_SAMPLES, qt, pt);
	save_data(qt, "Out/Pendulum/RefConv/pendulum_qt", N, split_order, stoch_order, seed);
	save_data(pt, "Out/Pendulum/RefConv/pendulum_pt", N, split_order, stoch_order, seed);
	save_data(sol.Brownian_process[0], "/gpfs/ashkerin/Conv_tests/Pendulum/RefProcess/BrownianProcess0", N, split_order, stoch_order, seed);
	save_data(sol.Brownian_process[1], "/gpfs/ashkerin/Conv_tests/Pendulum/RefProcess/BrownianProcess1", N, split_order, stoch_order, seed);
	save_data(sol.Brownian_process[2], "/gpfs/ashkerin/Conv_tests/Pendulum/RefProcess/BrownianProcess2", N, split_order, stoch_order, seed);

	return 0;
}