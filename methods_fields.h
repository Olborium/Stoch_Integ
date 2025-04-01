#ifndef METHODS_H
#define METHODS_H
#define _USE_MATH_DEFINES
#include <functional>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <map>
#include <stdexcept>
#include <fftw3.h>

extern double TEMP;             // temperature (in units of m^3/l, where l is the quartic coupling, m is the field mass)	
extern double SIZE;             // lattice size	(in units of inverse field mass)				
extern double DX;               // lattice step		
extern int N_x;                 // number of lattice points; x_{N_x} \equiv x_0; N_x must be power of 2.				
extern double DT;               // time step	

class BasicObjects {

    public:

        double* var_k;
        std::vector<double*> vars;                  // pointers to the variables [field, momentum]
        std::vector<fftw_plan> fft_f;               
        std::vector<fftw_plan> fft_i;
        std::vector<fftw_plan> fft_xk;

        BasicObjects();
        ~BasicObjects();
};

class InitialState : public BasicObjects {

    public:

        int sign;
        void prepare_state();

        InitialState(int sign);
};

class Evolve : public InitialState {

    public:
        
        double timespan;                            // simulation time
        double eta;                                 // damping coefficient
        int split_order;                            // order of splitting
        int stoch_order;                            // strong order

        Evolve(double eta, int split_order, int stoch_order, int sign); 

        void evolve(std::vector<double>&, double);

        double e_kin();                             // kinetic energy
        double e_pot();                             // potential energy
        void ps(std::vector<double>&, int);         // power spectrum
};

// save data to file
template <typename... Args>
void save_data(int prec, std::vector<double>& data, const std::string& name, Args&&... args)
{
	std::string long_name = name;
	((long_name += "_" + std::to_string(args)), ...);
	long_name += ".txt";	
	std::ofstream out(long_name, std::ios::app);
	if (!out.is_open()) {
		std::cerr << "Cannot open file!" << std::endl;
	} else {
		out << std::fixed << std::setprecision(prec);
        for(double& element : data){
			out << element << " ";
		}
		out << std::endl;
		out.close();
	}
};

#endif
