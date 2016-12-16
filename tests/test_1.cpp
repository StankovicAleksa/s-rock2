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
		const real x0=0.0;
 
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
					fy[i]=0.25*y[i]+0.5*sqrt(y[i]*y[i]+1);
	}

	real g_diag(const real x, const real yk, int k)
	{
			return sqrt((yk*yk+1)*0.5);
	}
  real init_condition(int k){
		return x0;
	}

	real rho(const real t, const real y[])
	{
			return abs(0.25+0.5*y[0]/sqrt(y[0]*y[0]+1));
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
  cout<<"\n----------1-dim test_1 ---------- "<<endl;
	simulation(argc,argv);
}

