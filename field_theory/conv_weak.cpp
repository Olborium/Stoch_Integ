//////////////////////////////
//  Weak convergence test   //
//////////////////////////////
//
// The script runs the real-time stochastic simulation of the scalar field with quartic self-interaction (stable potential).
// The field and momentum power spectra are measured at N_SAMPLE points during the run.
// The simulation time is TIME_SPAN * N_SAMPLE.
// The power spectra are saved to files.

#include "methods_fields.h"

int main(int argc, char* argv[])
{	
	SIZE = std::stod(argv[1]);		
	N_x = std::stoi(argv[2]);
	DT = std::stod(argv[3]);	
	int N_SAMPLE = std::stoi(argv[4]);
	double TIME_SPAN = std::stod(argv[5]);
	int split_order = std::stoi(argv[6]);
	int stoch_order = std::stoi(argv[7]);
	TEMP = std::stod(argv[8]);	
	double eta = std::stod(argv[9]);
	int sign = std::stoi(argv[10]);
	DX = SIZE/N_x;
///////////////////////////////////////////////////// 	
	std::vector<double> decay_times;
	std::vector<std::vector<double> > power_spectrum_phi(N_SAMPLE, std::vector<double>(N_x));
	std::vector<std::vector<double> > power_spectrum_chi(N_SAMPLE, std::vector<double>(N_x));	
	std::vector<double> phi_k2 (N_x);
	std::vector<double> chi_k2 (N_x);	
/////////////////////////////////////////////////////

//  initialize solver:
	Evolve sol(eta, split_order, stoch_order, sign);

//  sol.prepare_state() creates the thermal initial state, 
//	otherwise create the ground state:
	for(int j = 0; j < N_x; j++) {
		sol.vars[0][j] = 0.0;
		sol.vars[1][j] = 0.0;
	}
//  evolve and measure the power spectrum of the field and momentum:
	for (int n = 0; n < N_SAMPLE; n++) {	
		sol.evolve(decay_times, TIME_SPAN);
		sol.ps(phi_k2, 0);
		sol.ps(chi_k2, 1);
		for (int j = 0; j < N_x; j++) {
			power_spectrum_phi[n][j] = phi_k2[j];
			power_spectrum_chi[n][j] = chi_k2[j];
		}
	}
//  save the data:	
	std::random_device rd;
	int seed = rd();
	for (auto& time_slice : power_spectrum_phi) {
		save_data(5, time_slice, "/home/olborium/scratch/power_spectrum_phi/ps_", eta, TEMP, DT, split_order, stoch_order, seed);
	}	
	for (auto& time_slice : power_spectrum_chi) {
		save_data(5, time_slice, "/home/olborium/scratch/power_spectrum_chi/ps_", eta, TEMP, DT, split_order, stoch_order, seed);
	}	
	return 0;
}
