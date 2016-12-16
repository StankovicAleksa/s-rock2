#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <iostream>
#include <cmath>
#include <set>
namespace parameters{
    typedef double real;
    
    extern real uround;
    
    // Number of equations
    extern int neqn;

		// Universal parameters
    extern real tbeg; //starting time
    extern real tend; //end time
		
		//-----------------------------------------------------------
		// Equation of the problem:
		//
		// dX=f(X)dt + g(X)dW
		// X(0) = init_condition
		//
		//-----------------------------------------------------------
		void f(const real x, const real y[], real fy[]);
		// stochastic term
		real g_diag(const real x, const real yk, int k);
		// initial condition
		real init_condition(int k);

		//-----------------------------------------------------------
		//    The subroutine RHO gives an estimation of the spectral 
		//    radius of the Jacobian matrix of the problem. This
		//    is a bound for the whole interval and thus RHO is called
		//    once.
		//-----------------------------------------------------------
		real rho(const real t, const real y[]);
	
		// Additional functions 	
		void init_parameters();
    double get_cpu_time();
		

		extern long long f_evals;
		extern long long g_evals;
}

#endif
