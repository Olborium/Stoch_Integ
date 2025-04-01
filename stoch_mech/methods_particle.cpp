#include "methods_particle.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
							Constants
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Global constants:

lfloat TEMP; // temperature
lfloat DT;   // time step

// Splitting coefficients:

std::vector<std::vector<lfloat>> A = {
{1.0Q},
{0.5Q, 0.5Q},
{0.29166666666666666666666666666667Q, 0.75Q, -0.04166666666666666666666666666667Q},
{0.1344961992774310892Q, -0.2248198030794208058Q, 0.7563200005156682911Q, 0.3340036032863214255Q}
};
std::vector<std::vector<lfloat>> B = {
{1.0Q},
{1.0Q, 0.0Q},
{0.666666666666666666666666666666667Q, -0.666666666666666666666666666666667Q, 1.0Q},
{0.5153528374311229364Q, -0.085782019412973646Q, 0.4415830236164665242Q, 0.1288461583653841854Q}
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
							Noise generator
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct S { lfloat d1; lfloat d2; };

struct Q { lfloat d1; lfloat d2; lfloat d3; };

void lprint(lfloat a) {
    char buffer[128];
    quadmath_snprintf(buffer, sizeof(buffer), "%.36Qg", a);
    std::cout << buffer << std::endl;	
}

void generate_noise( std::vector<Q>& q, bool strong_conv,
					 std::vector<std::vector<lfloat>>& Brownian_process, lfloat dt, int stoch_order)
{
	int num_of_steps = q.size();
	int ord_st[3] = {1, stoch_order > 1 ? 1 : 0, stoch_order > 2 ? 1 : 0};
	if ( strong_conv==false ) { 		
// Calculate independent noise variables \zeta_{1,2,3}: for weak convergence tests and practical implementation
		std::random_device rd;
		std::mt19937 gen(rd());
		std::normal_distribution<long double> distr(0.0L, 1.0L);		
		for (int i = 0; i < num_of_steps; i++ ) {
			lfloat n1 = static_cast<lfloat>(distr(gen));
			lfloat n2 = static_cast<lfloat>(distr(gen));
			lfloat n3 = static_cast<lfloat>(distr(gen));
			q[i].d1 = ord_st[0]*sqrtq(1.0/DT)*n1;
			q[i].d2 = ord_st[1]*sqrtq(DT/12.0)*n2;
			q[i].d3 = ord_st[2]*DT*sqrtq(DT/720.0)*n3;	
		}				
	} else { 
// Calculate noise variables \zeta_{1,2,3} from the underlying Brownian process: for strong convergence tests
		int points_per_step = (Brownian_process[0].size()-1)/num_of_steps;
		for (int i = 0; i < num_of_steps; i++) {
			q[i].d1 = 1.0/DT*(Brownian_process[0][(i+1)*points_per_step] - Brownian_process[0][i*points_per_step]);
			q[i].d2 = 1.0/DT*(Brownian_process[1][(i+1)*points_per_step] - Brownian_process[1][i*points_per_step])
					+0.5*(1.0-dt/DT)*(Brownian_process[0][(i+1)*points_per_step] - Brownian_process[0][i*points_per_step]);
			q[i].d3 = 1.0/DT*(Brownian_process[2][(i+1)*points_per_step] - Brownian_process[2][i*points_per_step])
					+0.5*(1.0-dt/DT)*(Brownian_process[1][(i+1)*points_per_step] - Brownian_process[1][i*points_per_step])	
					+ (DT/12.0-dt/4.0+dt*dt/6.0/DT) * (Brownian_process[0][(i+1)*points_per_step] - Brownian_process[0][i*points_per_step]);
			for (int j = 0; j < points_per_step; j++) {
				lfloat tj = j*dt;
				lfloat dw0 = Brownian_process[0][i*points_per_step+j+1] - Brownian_process[0][i*points_per_step+j];
				lfloat dw1 = Brownian_process[1][i*points_per_step+j+1] - Brownian_process[1][i*points_per_step+j];				
				q[i].d2 += - tj/DT * dw0;
				q[i].d3 += - tj/DT * dw1 + (-tj/2.0+tj*tj/2.0/DT+tj*dt/2.0/DT) * dw0;
			}
			q[i].d1 = ord_st[0]*q[i].d1;
			q[i].d2 = ord_st[1]*q[i].d2;
			q[i].d3 = ord_st[2]*q[i].d3;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
							Evolution
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Integrator with PQ splitting
void evol_PQ( lfloat* vars, bool strong_conv, std::vector<std::vector<lfloat>>& Brownian_process, lfloat dt,
			  int num_of_steps, lfloat eta, int split_order, int stoch_order, int num_of_samples, 
			  std::vector<lfloat>& qt, std::vector<lfloat>& pt )
{
	S s;
	std::vector<Q> q(num_of_steps);
	lfloat sigma = sqrtq(2.0*eta*TEMP);
	int w = split_order - 1;
	
	generate_noise( q, strong_conv, Brownian_process, dt, stoch_order );

	for (int i = 0; i < num_of_steps; i++) {
	// Calculate stochastic forces	
		s.d1 = q[i].d2 - eta*q[i].d3;
		s.d2 = q[i].d1 - eta*q[i].d2 + (eta*eta-1.0-3.0*vars[0]*vars[0])*q[i].d3;    // anharmonic oscillator
//		s.d2 = q[i].d1 - eta*q[i].d2 + (eta*eta-cosq(vars[0]))*q[i].d3;    			 // pendulum
		for (int n = 0; n < split_order; n++) {
	// Solve 1st equation (P)
			lfloat F = -vars[0]-vars[0]*vars[0]*vars[0]+sigma*s.d2;                  // anharmonic oscillator
//			lfloat F = -sinq(vars[0])+sigma*s.d2;                                    // pendulum
			if(eta==0.0L) {
				vars[1] = vars[1] + A[w][n]*DT*F;
			} else {
				vars[1] = F/eta + expq(-eta*A[w][n]*DT)*(vars[1]-F/eta);
			}
	// Solve 2nd equation (Q)
			vars[0] = vars[0] + B[w][n]*DT*(vars[1]+sigma*s.d1);				
		}
	// write to the output
		if( (i+1)%(num_of_steps/num_of_samples) == 0 ) {
			qt.push_back(vars[0]);
			pt.push_back(vars[1]);
		}		
	}
}

// Integrator with LN splitting
void evol_LN( lfloat* vars, bool strong_conv, std::vector<std::vector<lfloat>>& Brownian_process, lfloat dt,
			  int num_of_steps, lfloat eta, int split_order, int stoch_order, int num_of_samples,
			  std::vector<lfloat>& qt, std::vector<lfloat>& pt )
{
	S s;
	std::vector<Q> q(num_of_steps);
	lfloat sigma = sqrtq(2.0*eta*TEMP);
	int w = split_order - 1;		
	bool overdamped;
	if (1.0Q > eta/2.0Q) { overdamped = false; } else { overdamped = true; }	

	generate_noise( q, strong_conv, Brownian_process, dt, stoch_order );

	for (int i = 0; i < num_of_steps; i++) {
	// Calculate stochastic forces
		s.d1 = q[i].d2 - eta*q[i].d3;
 		s.d2 = q[i].d1 - eta*q[i].d2 + (eta*eta-1.0Q-3.0Q*vars[0]*vars[0])*q[i].d3; // anharmonic oscillator
		for (int n = 0; n < split_order; n++) {
	// Kick (N)
			vars[1] = vars[1] - A[w][n]*DT*vars[0]*vars[0]*vars[0];                 // anharmonic oscillator
	// Drift (L)			
			lfloat Oe = sqrtq(fabsq(1.0Q-eta*eta/4.0Q));
			lfloat a0 = vars[0] - sigma*(eta*s.d1+s.d2);
			lfloat b0 = (vars[1]+eta*vars[0]/2.0Q)/Oe + sigma*(s.d1-eta*(eta*s.d1+s.d2)/2.0Q)/Oe;
			lfloat c0 = sigma*(eta*s.d1+s.d2);
			lfloat a1 = vars[1] + sigma*s.d1;
			lfloat b1 = -(vars[0]+eta*vars[1]/2.0Q)/Oe + sigma*(eta*s.d1+2.0Q*s.d2)/2.0Q/Oe;
			lfloat c1 = -sigma*s.d1;
			if (overdamped == false) {
				vars[0] = (a0*cosq(B[w][n]*DT*Oe) + b0*sinq(B[w][n]*DT*Oe))*expq(-eta*B[w][n]*DT/2.0Q) + c0;
				vars[1] = (a1*cosq(B[w][n]*DT*Oe) + b1*sinq(B[w][n]*DT*Oe))*expq(-eta*B[w][n]*DT/2.0Q) + c1;
			} else {
				vars[0] = (a0*coshq(B[w][n]*DT*Oe) + b0*sinhq(B[w][n]*DT*Oe))*expq(-eta*B[w][n]*DT/2.0Q) + c0;
				vars[1] = (a1*coshq(B[w][n]*DT*Oe) + b1*sinhq(B[w][n]*DT*Oe))*expq(-eta*B[w][n]*DT/2.0Q) + c1;				
			}					
		}
	// write to the output
		if( (i+1)%(num_of_steps/num_of_samples) == 0 ) {
			qt.push_back(vars[0]);
			pt.push_back(vars[1]);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
							Methods of classes
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BasicObjects::BasicObjects( lfloat q0, lfloat p0 ) : q0(q0), p0(p0) {}

BasicObjects::~BasicObjects() = default;	  

Evolve::Evolve( lfloat q0, lfloat p0, lfloat eta, lfloat dt, const std::string& INTEG, bool strong_conv, int split_order, int stoch_order)
	try : eta(eta), dt(dt), INTEG(INTEG), split_order(split_order), stoch_order(stoch_order), strong_conv(strong_conv),
		  functionMap({
			{"LN", evol_LN},
			{"PQ", evol_PQ}
		  }), BasicObjects(q0, p0) {
			if (functionMap.count(INTEG) == 0) {
				throw std::invalid_argument("Invalid splitting type!");
			}
		  }
	catch (const std::invalid_argument& e) {
		std::cerr << e.what() << std::endl;
	} 	

Evolve::~Evolve() = default;

void Evolve::evolve(int num_of_steps, int num_of_samples, std::vector<lfloat>& qt, std::vector<lfloat>& pt) {
	return Evolve::functionMap[INTEG]( BasicObjects::vars, Evolve::strong_conv, Evolve::Brownian_process,
	 Evolve::dt, num_of_steps, Evolve::eta, Evolve::split_order, Evolve::stoch_order, num_of_samples, qt, pt );
}

void Evolve::generate_Brownian_process( int num_of_steps, lfloat dt, int seed ) 
{
	long double ddt = static_cast<long double>(dt);
	std::mt19937 gen(seed);
	Evolve::Brownian_process.push_back({0.0Q});
	Evolve::Brownian_process.push_back({0.0Q});
	Evolve::Brownian_process.push_back({0.0Q});
	std::normal_distribution distr0(0.0L,std::sqrt(ddt));
	std::normal_distribution distr1(0.0L,std::sqrt(ddt*ddt*ddt/12.0));
	std::normal_distribution distr2(0.0L,std::sqrt(ddt*ddt*ddt*ddt*ddt/720.0));
	for(int i = 1; i < num_of_steps+1; i++) {		
		lfloat b0 = Evolve::Brownian_process[0].back();
		lfloat b1 = Evolve::Brownian_process[1].back();
		lfloat b2 = Evolve::Brownian_process[2].back();
		lfloat d0 = static_cast<lfloat>(distr0(gen));
		lfloat d1 = static_cast<lfloat>(distr1(gen));
		lfloat d2 = static_cast<lfloat>(distr2(gen));
		Evolve::Brownian_process[0].push_back( b0 + d0 );
		Evolve::Brownian_process[1].push_back( b1 + d1 );
		Evolve::Brownian_process[2].push_back( b2 + d2 );
	}
}

void Evolve::load_Brownian_process( const std::string& path, int num_of_steps, lfloat dt, int seed )
{
	Evolve::Brownian_process.push_back({});
	Evolve::Brownian_process.push_back({});
	Evolve::Brownian_process.push_back({});
	for( int i = 0; i < 3; i++ ) {
		std::string name = path + "/BrownianProcess"+std::to_string(i)+"_"
		+std::to_string(num_of_steps)+"_2_1_"
		+std::to_string(seed)+".txt";
		std::ifstream file(name);
		if (!file.is_open()) {
			std::cerr << "Cannot open file!" << std::endl;
		} else {
			std::string line;
			std::getline(file, line);
			std::stringstream ss(line);
			std::string numStr;
			while (ss >> numStr) {
				const char* c_numStr = numStr.c_str();
				lfloat value = strtoflt128(c_numStr, NULL);
				Evolve::Brownian_process[i].push_back(value);
			}
			file.close();
		}
	}
}  
