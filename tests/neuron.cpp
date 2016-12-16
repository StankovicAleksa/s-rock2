//* * * * * * * * * * * * * * * * * * * * * * * * *
//   Driver for ROCK4 (or ROCK2) at Brusselator-2dim problem
//* * * * * * * * * * * * * * * * * * * * * * * * *
//
//   This driver shows how to use S-ROCK2. It solves a
//   system of SDEs resulting from the 1-dimensional space 
//   discretization of the Neuron cable equation.
//    
//   du= nu*u_xx dt - beta*u dt + sigma*(u+v0)/dx dW(t)
//
//   with initial conditions 
//
//   u(x,0)=-70+20cos(5pi*x)(1-x)
//
//   and homogeneous Neumann boundary conditions
//
//   u_x(0,t)=u_x(1,t)=0.
//
//   We discretize the space variables with
//   x_i=i/(N-1) for i=0,1,...,N-1,
//   The spectral radius of the Jacobian 
//   can be estimated with the Gershgorin theorem   
//   by 4*(N-1)^2. Thus we provide an external
//   function RHO, giving the spectral radius 
//   of the Jacobian matrix.
//   As output point we choose t_out=1.
//
//-------------------------------------------------------

#include <iostream>
#include <limits>
#include <random>
#include "../integrators/parameters.h"
#include "simulation.h"
#include <cmath>
using namespace std;
// set parameters
namespace parameters{
    typedef double real;
    
    real uround;
    // Number of equations
    int neqn=128;

		// Universal constants
    real tbeg=0.0; //starting time
    real tend=1.0; //end time
    
    // Problem constants
    real nu=0.01;
    real beta=1.;
    real sigmaeq=4.*1e-3;
    real V0=10.;
    real dx;
    real sqrdx;
    const real pi=4.*atan(1.);

 
    void init_parameters()
    {
        uround=std::numeric_limits<real>::epsilon();
        dx= 1./(neqn-1.);
        sqrdx=sqrt(dx);
    }
    
    double get_cpu_time(){
        return (double)clock() / CLOCKS_PER_SEC;
    };

	// Computes the deterministic part of the righthand side of the SDE
	void f(const real t, const real y[], real fy[])
	{
			// Computing diffusion with Neumann bnd conditions
			fy[0]=nu*2.*(y[1]-y[0])*(neqn-1.)*(neqn-1.)-beta*y[0];
			fy[neqn-1]=nu*2.*(y[neqn-2]-y[neqn-1])*(neqn-1.)*(neqn-1.)-beta*y[neqn-1];
			for (int i=1;i<neqn-1;i++)
			{
					fy[i]=nu*(y[i-1]-2.*y[i]+y[i+1])*(neqn-1.)*(neqn-1.)-beta*y[i];
					if(abs(i/(neqn-1.)-0.5)<0.1)
							fy[i] += 5.*exp(1.-1e-2/(1e-2-(i/(neqn-1.)-0-5)*(i/(neqn-1.)-0.5)));
			}
	}

	real g_diag(const real x, const real yk, int k)
	{
			return sigmaeq*(yk+V0)/sqrdx;
	}
  real init_condition(int k){
    real x=((real)k)/(neqn-1.);
		return -70.+20.*cos(15.*pi*x)*(1.-x);
	}

	real rho(const real t, const real y[])
	{
			return nu*4.*(neqn-1.)*(neqn-1.)+beta;
	}
	real get_integral(real y[])
	{
			real dx = 1./(neqn-1.);
			real integ;
			
			integ = 0.5*dx*y[0]+0.5*dx*y[neqn-1];
			for(int i=1;i<neqn-1;i++)
					integ += dx*y[i];
			
			return integ;
	}
}

// run equation
int main(int argc, char *argv[]){
  cout<<"\n----------1-dim Neuron---------- "<<endl;
	simulation(argc,argv);
}

