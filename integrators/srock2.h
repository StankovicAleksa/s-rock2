#include <random>

#include "parameters.h"
#include <algorithm>
using namespace std;
using namespace parameters;

namespace rock2{
		struct stoch_sample{
			double t;
			double val;
			stoch_sample(double _t,double _val){t=_t;val=_val;}
		};	
    extern int ms[46];
    extern real fp1[46];
    extern real fp2[46];
    extern real recf[4476];
    extern real recalph[46];
    extern real recf2[184];
    extern std::mt19937 generator; // randomization generator
    extern uniform_int_distribution<> distribution;
    void mdegr(int& mdeg, int mp[]);

		// 
    void rock2( real& t, const real tend, real& h, real*& y,
              real*& yn, real*& yjm1, real*& yjm2, real*& fn, real*& stowork,
              const real atol[], const real rtol[], 
              int iwork[], int& idid);

    void rtstep_no_noise( const real t, const real h, real*& y,
                real*& yn, real*& fn, real*& yjm1, real*& yjm2,
                int mdeg, int mp[], real& err, const real atol[],
                const real rtol[], const int ntol);
    
    void rtstep_diag_noise(const real t, const real h, real*& y,
                real*& yn, real*& fn, real*& yjm1, real*& yjm2, real*& stowork,
                int mdeg, int mp[], real& det_err,real& stoch_err, const real atol[],
                const real rtol[], const int ntol);
    
    void rtstep_gen_noise( const real t, const real h, real*& y,
                real*& yn, real*& fn, real*& yjm1, real*& yjm2,
                int mdeg, int mp[], real& err, const real atol[],
                const real rtol[], const int ntol);
    
    void rtstep_no_noise_rkc(const real t, const real h, real*& y,
            real*& yn, real*& fn, real*& yjm1, real*& yjm2,
            int mdeg, real& err, const real atol[],
            const real rtol[], const int ntol);

    void check_initialisation(const real& t, const real tend, real& h,
               real atol[], real rtol[], 
               int iwork[], int& idid);
    
    real rocktrho(const real t,
                     real *yn, real *fn, real *z, real *fz, real *pz,
                     int& idid, int iwork[]);


		/*** Required for timestepping ****/
		struct sp_sample{
			real val;
			real t;
			sp_sample(real _t, real _val){ t=_t; val=_val;}
			sp_sample(){ t=0; val=0;} // default constructor
		};
		bool operator < (sp_sample s1, sp_sample s2);

		real get_Jqq(const real &xi,const real &h);
  	void construct_xi(const real &t,const real &h, const int& mdeg,real *xi);  
		void construct_chi(const real &t, const real &h, const int& mdeg, real *chi);  // Added by Aleksa
  	void construct_Jqr(const int &mdeg,const real &h, const real*xi, const real *chi,real **&J_qr); // Added by Aleksa
		extern std::set<sp_sample> *chi_t;
		extern std::set<sp_sample> *xi_t;
}
