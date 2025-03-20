## Stoch_Integ

The code and testing routines for the paper 2503.xxxxx.

# Strong convergence tests

- Main routines:
  
  methods_particle.h, methods_particle.cpp

- Generate 10 realizations of white noise and 'true' solution:
  
  p_conv_strong_ref_proc.sbatch running p_conv_strong_ref_proc.cpp

- Generate approximate solutions using the reference white noise and high-order stochastic schemes; measure the convergence rate to the 'true' solution:
  
  p_ref_conv.sh running p_ref_conv.cpp, p_ref_conv.py

- Generate approximate solutions using the reference white noise and high-order stochastic schemes; measure the rate of self-convergence:
  
  p_self_conv.sh running p_self_conv.cpp, p_self_conv.py

# Application to field theory

- Main routines:
  
  methods_fields.h, methods_fields.cpp

- Run a series of real-time stochastic simulations of the real scalar field with quartic self-interaction (stable potential) in 1+1 dimensions;
  
  measure the thermodynamic observables and compare them with the theory predictions; measure the rate of weak convergence with respect to these observables:
  
  conv_weak.sbatch running conv_weak.cpp, conv_weak_temp.py, conv_weak_power_spectrum.py

- Run a series of real-time stochastic simulations of the real scalar field with quartic self-interaction (unstable potential) in 1+1 dimensions;
