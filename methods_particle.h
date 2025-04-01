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
#include <quadmath.h>           // the library implementing type float_128

using lfloat = __float128;      // this type guarantees 128 bit (~33 digits) precision

extern lfloat TEMP;             // temperature					
extern lfloat DT;               // time step	

void lprint(lfloat);

class BasicObjects {

    public:

        lfloat q0;              // initial coordinate
        lfloat p0;              // initial momentum
        lfloat vars[2] = {q0, p0}; 

        BasicObjects(lfloat, lfloat);
        ~BasicObjects();
};

class Evolve : public BasicObjects {

    private:

        std::string INTEG;     

        std::map<std::string, std::function<void( lfloat*, bool, std::vector<std::vector<lfloat>>&, lfloat, int, lfloat, 
                                                  int, int, int, std::vector<lfloat>&, std::vector<lfloat>& ) > > functionMap;
    
    public:

		lfloat eta;             // damping coefficient
        lfloat dt;              // time step of the reference Brownian process
        int split_order;        // order of splitting of the symplectic scheme
        int stoch_order;        // order of strong convergence
        bool strong_conv;       // strong convergence test: yes or no

        std::vector<std::vector<lfloat>> Brownian_process;
        void generate_Brownian_process(int, lfloat, int);
        void load_Brownian_process(const std::string&, int, lfloat, int);
        
        Evolve(lfloat, lfloat, lfloat, lfloat, const std::string&, bool, int, int);      
        ~Evolve();

        void evolve(int, int, std::vector<lfloat>&, std::vector<lfloat>&);           
};

// Save data to file
template <typename... Args>
void save_data(std::vector<lfloat>& data, const std::string& name, Args&&... args)
{
	std::string long_name = name;
	((long_name += "_" + std::to_string(args)), ...);
	long_name += ".txt";	
	std::ofstream out(long_name, std::ios::app);
	if (!out.is_open()) {
		std::cerr << "Cannot open file!" << std::endl;
	} else {
        for(lfloat& element : data){
			char buffer[128];
            quadmath_snprintf(buffer, sizeof(buffer), "%.36Qg", element);
            out << buffer << " ";
		}		
		out << std::endl;
		out.close();
	}
};

#endif