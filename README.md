## Stoch_Integ

The code and testing routines for the paper 2503.xxxxx.

# Strong convergence tests

- Main routines:
  
    _methods_particle.h, methods_particle.cpp_

- Generate 10 realizations of white noise and 'true' solution:
  
    _p_conv_strong_ref_proc.sbatch -> p_conv_strong_ref_proc.cpp_

- Generate approximate solutions using the reference white noise and high-order stochastic schemes; measure the convergence rate to the 'true' solution:
  
    _p_ref_conv.sh -> p_ref_conv.cpp, p_ref_conv.py_

- Generate approximate solutions using the reference white noise and high-order stochastic schemes; measure the rate of self-convergence:
  
    _p_self_conv.sh -> p_self_conv.cpp, p_self_conv.py_

# Application to field theory

- Main routines:
  
    _methods_fields.h, methods_fields.cpp_

- Run a series of real-time stochastic simulations of the real scalar field with quartic self-interaction (stable potential) in 1+1 dimensions;
  measure the thermodynamic observables and compare them with the theory predictions; measure the rate of weak convergence with respect to these observables:
  
    _conv_weak.sbatch -> conv_weak.cpp, conv_weak_temp.py, conv_weak_power_spectrum.py_

- Run a series of real-time stochastic simulations of the real scalar field with quartic self-interaction (unstable potential) in 1+1 dimensions;
  measure the decay times from the initial state of thermal equilibrium around the false vacuum; measure the decay rate from the solve of the survival probability:

  _decay_rate.sbatch -> decay_rate.cpp_
  _decay_rate.py_
