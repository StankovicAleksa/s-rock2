#include <iostream>
#include <limits>
#include <random>
#include "../integrators/parameters.h"
#include "simulation.h"
#include <cmath>
// set parameters
using namespace std;
namespace parameters{
    typedef double real;
    
    real uround;
    // Number of equations
    int neqn=1;

		// Universal constants
    real tbeg=0.0; //starting time
    real tend=1.0; //end time
    
    // Problem constants
		const real lambda=1.0;
		const real mu=0.01;
    const real pi=4.*atan(1.);
		const real x0=1.0;
 
    void init_parameters()
    {
        uround=std::numeric_limits<real>::epsilon();
    }
    
    double get_cpu_time(){
        return (double)clock() / CLOCKS_PER_SEC;
    };

	// Computes the deterministic part of the righthand side of the SDE
	void f(const real t, const real y[], real fy[])
	{
			for (int i=0;i<neqn;i++)
					fy[i]=lambda*y[i];
	}

	real g_diag(const real x, const real yk, int k)
	{
			return mu*yk;
	}
  real init_condition(int k){
		return x0;
	}

	real rho(const real t, const real y[])
	{
			return abs(lambda);
	}
	//-----------------------------------------------------------------
	// Calculate exact solution
	//-----------------------------------------------------------------
	real get_integral(real y[])
	{
			real integ=0.0;
			for(int i=0;i<neqn;i++)
					integ += y[i];
			
			return integ;
	}
}

// run equation
int main(int argc, char *argv[]){
  cout<<"\n----------1-dim Black-Scholes ---------- "<<endl;
	simulation(argc,argv);
}

