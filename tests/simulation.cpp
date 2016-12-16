#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <time.h>
#include <sys/time.h>

#include "simulation.h"
#include "../integrators/srock2.h"
using namespace rock2;
using namespace std;
namespace parameters{
	long long f_evals,g_evals;
}
void output_sol(const real x[], const real t,
                const int nsol=1, string solfile="sol.m");
void output_mesh(string solfile="data/mesh.m");
void output_means(real mean, string mtfile="mtdata.m");
void output_evals_f(int evals, string evals_f_file);
void output_evals_g(int evals, string evals_g_file);

void output_statistics(string stat_file,real precision, real mean_int, int f_evals, int g_evals);
void square(real x[], real y[]);
void exp(real x[], real y[]);
void sin(real x[], real y[]);
void arcsinh_2(real x[], real y[]);


int simulation(int argc, char *argv[])
{
    //-------------------------------------------------------
    //    Initialize iwork: 
    //    If iwork(0)=1  RHO returns an upper bound for the spectral radius
    //    else           it is computed internally by rock2.
    //    If iwork(1)=1  The Jacobian is constant (RHO is called once).
    //    If iwork(2)=0  Return a solution at tend.
    //    else           return a solution after every step.
    //    If iwork(3)=0  Atol and rtol are scalars.
    //    If iwork(12)=1 Use adaptive time step.
    //    If      iwork(13)=0  No noise (deterministic problem),
    //    else if iwork(13)=1  diagonal noise,
    //    else if iwork(13)=2  general noise.
    //    else if iwork(13)=3  no noise but with RKC instead of ROCK2
    //    If iwork[14]=1       Verbose
    //-------------------------------------------------------
    int iwork[17],idid;
    iwork[0]=1;
    iwork[1]=0;
    iwork[2]=1;
    iwork[3]=0;
    iwork[12]=0;
    iwork[13]=1;
    iwork[14]=0;
		f_evals=g_evals=0;
    real dt;//time step
    int mt; //monte carlo iterations
    int p; //processor
    
    if(argc>1){//reading from console
        istringstream hbuf(argv[1]);
        istringstream mtbuf(argv[2]);
        istringstream pbuf(argv[3]);
        hbuf>>dt;
        mtbuf>>mt;
        pbuf>>p; }
    else{//or standard parameters
        dt=1e-2;
        mt = 1;
        p=1;
    }
		mt=10000;
    const real cdt=dt; //needed for restarting iterations

    //----- required tolerance (needed for adaptive time step) -----
    real rtol[1]={1.0e-2};
    real atol[1]={rtol[0]};

    //output file
    ostringstream f_sstr;
		// If adaptive timestep is enabled we store results in files
		// named with respect to tolerances
		if ( iwork[12] == 1 ){ 
    	f_sstr<<scientific<<setprecision(0)<<atol[0];
		}
		else{
    	f_sstr<<scientific<<setprecision(0)<<dt;
		}
    string dt_str = f_sstr.str();
    ostringstream p_sstr;
    p_sstr<<p;
    string p_str = p_sstr.str();
    string solfile="data/sol.m";
    string mtfile_id="data/mtdata_id_"+dt_str+"_"+p_str+".txt";
    string mtfile_sq="data/mtdata_sq_"+dt_str+"_"+p_str+".txt";
    string mtfile_sin="data/mtdata_sin_"+dt_str+"_"+p_str+".txt";
    string mtfile_exp="data/mtdata_exp_"+dt_str+"_"+p_str+".txt";
    string mtfile_arcsinh_2="data/mtdata_arcsinh_2_"+dt_str+"_"+p_str+".txt";
   	
		string evalsfile_f ="data/evals_f_"+dt_str+"_"+p_str+".txt" ;
		string evalsfile_g ="data/evals_g_"+dt_str+"_"+p_str+".txt" ;

    init_parameters();
    
    // all the variables will be contained in the work vectors
    real work[6*neqn];
    real* stowork = new real[7*neqn];
		// Allocate distribution variables
		chi_t=new set<rock2::sp_sample>[neqn]; 		
		xi_t =new set<sp_sample>[neqn]; 		
		// pointer to the solution and working vectors
    real *y=work;
    real *yn=work+neqn;
    real *yjm1=work+2*neqn;
    real *yjm2=work+3*neqn;
    real *fn=work+4*neqn;
    real *ymean=work+5*neqn;
    for(int i=0;i<neqn;i++)
        ymean[i]=0;
    

    rock2::check_initialisation(tbeg, tend, dt, atol, rtol, iwork, idid);

    //----- integration -----
    cout<<"\nProblem Configuration:"<<endl;
    cout<<"Mesh dimension: "<<neqn;
    cout<<"Starting step-size dt= "<<dt<<endl;
    cout<<"Time-step adaptivity enabled: "<<(iwork[12]==1 ? "True":"False")<<endl;
    cout<<"Spectral radius computed "<<(iwork[0]==0 ? "internally":"by rho")<<endl;
    cout<<"Spectral radius is constant: "<<(iwork[1]==0 ? "False":"True")<<endl;
    cout<<"Initial time t_0= "<<tbeg<<endl;
    cout<<"End time t_e= "<<tend<<endl;
    cout<<"ATol = RTol = "<<rtol[0]<<endl;
    cout<<"Solution is tabulated in file sol.out"<<endl;
    
    real cpustart,cpufinish;
    int nsol;
    const int foutput=3000; //frequency output
        
    //monte carlo variables
    real mean_int, mean_int_sq, mean_int_sin, mean_int_exp, mean_int_arcsinh_2;
    mean_int=mean_int_sq=mean_int_sin=mean_int_exp=mean_int_arcsinh_2=0.;
    
    
    
    cpustart=get_cpu_time(); //measuring time makes sense if we are not making
                             //outputs in the loop, i.e. iwork[2]=0
    
    
    // begin monte carlo iterations
    for(int mtiter=1;mtiter<=mt;mtiter++)
    {
        //not really needed
        int iter=0;
        nsol=1;

				// added
    		init_parameters();
				y=work;
				yn=work+neqn;
				yjm1=work+2*neqn;
				yjm2=work+3*neqn;
				fn=work+4*neqn;
				ymean=work+5*neqn;
				for (int i=0;i<6*neqn;i++) work[i]=0;
				for (int i=0;i<7*neqn;i++) stowork[i]=0;
    		for(int i=0;i<neqn;i++)
       		 ymean[i]=0;
    		rock2::check_initialisation(tbeg, tend, dt, atol, rtol, iwork, idid);
        
        if(mtiter%foutput==0)
            cout<<"Monte Carlo iteration "<<mtiter<<endl;
        
        //reinitialize statistics vector
        for(int i=4;i<12;i++)
            iwork[i]=0;
       
        //reinitialize y
        real t=tbeg;
        dt = cdt;
        const real pi=4.*atan(1.);
        for (int j=0;j<neqn;j++)
						y[j]=init_condition(j);

       	// delete history of sampling
       	for ( int i=0;i<neqn;i++){
					xi_t[i].clear(),chi_t[i].clear();
					xi_t[i].insert(sp_sample(tbeg,0));
					chi_t[i].insert(sp_sample(tbeg,0));  // set value at t=t_0
				}	
				
				if(mt==1)
            output_sol(y,t,nsol,solfile); //output initial value

        while(t<tend)
        {
            iter++;
            rock2::rock2(t,tend,dt,y,
                 yn,yjm1,yjm2,fn,stowork, //working vectors
                 atol,rtol,
                 iwork,idid);
            if(iter%foutput==0 && t<tend)
                output_sol(y,t,++nsol,solfile); //output current solution
        }
        mean_int += get_integral(y)/mt;
        sin(y,yn);
        mean_int_sin+=get_integral(yn)/mt;
        square(y,yn);
        mean_int_sq+=get_integral(yn)/mt;
        exp(y,yn);
        mean_int_exp+=get_integral(yn)/mt;
				arcsinh_2(y,yn);
				mean_int_arcsinh_2+=get_integral(yn);
        for(int i=0;i<neqn;i++)
            ymean[i] += y[i]/mt;
    }
		f_evals/=mt;
		g_evals/=mt;
		mean_int_arcsinh_2/=mt;
    cpufinish=get_cpu_time();
    
     //----- print data -----
    if(mt==1)
    {
        output_sol(ymean,tend,++nsol,solfile);
        output_mesh();
    }
		//output_means(mean_int_arcsinh_2,mtfile_arcsinh_2);
		if ( iwork[12] == 0 )
			output_statistics("data/stats.txt",dt,mean_int,f_evals,g_evals);
		else
			output_statistics("data/stats.txt",atol[0],mean_int,f_evals,g_evals);
		/*
    output_means(mean_int,mtfile_id);
    output_means(mean_int_sin,mtfile_sin);
    output_means(mean_int_sq,mtfile_sq);
    output_means(mean_int_exp,mtfile_exp);
   
	 	output_evals_f(f_evals,evalsfile_f);	
	 	output_evals_g(g_evals,evalsfile_g);	
    //------ save elapsed time
    string timefile = "data/time_"+p_str+".txt";
    ofstream cputimesave(timefile.c_str(),std::ofstream::out|std::ofstream::app);
    cputimesave<<cpufinish-cpustart<<endl;
    cputimesave.close();
		*/
    //----- print statistics -----
    cout<<"\nStatistics:"<<endl;
    cout<<"The value of IDID is "<<idid<<endl;
    cout<<"Max estimation of the spectral radius = "<<iwork[10]<<endl;
    cout<<"Min estimation of the spectral radius = "<<iwork[11]<<endl;
    cout<<"Max number of stages used = "<<iwork[9]<<endl;
    cout<<"Number of f evaluations for the spectral radius = "<<iwork[8]<<endl;
    cout<<"Number of f evaluations for integration ="<<f_evals<<endl;
    cout<<"Number of g evaluations for integration ="<<g_evals<<endl;
    cout<<"Steps: "<<iwork[5]<<endl;
    cout<<"Accepted steps: "<<iwork[6]<<endl;
    cout<<"Rejected steps: "<<iwork[7]<<endl; 
    cout<<"Number of Monte Carlo samples: "<<mt<<endl;
    cout<<"Elapsed time: "<<cpufinish-cpustart<<" secs"<<endl;
    delete[] stowork;
    
    //-------------------------------------------------------
    //    End of main program
    //-------------------------------------------------------
}


      
void output_sol(const real x[], const real t, 
                const int nsol, string solfile)
{
		// Quick fix, do not output sol( files grow too big)

    ofstream file;
    if(nsol==1)
        file.open(solfile.c_str(),std::ofstream::out);
    else
        file.open(solfile.c_str(),std::ofstream::out|std::ofstream::app);
    file.precision(16);
    
    // writing solution in matlab syntax
    file<<"t("<<nsol<<")="<<t<<";"<<endl;
    file<<"y(:,"<<nsol<<")=[";
    for(int j=0;j<neqn;j++)
    {
        file<<x[j];
        if(j<neqn-1) file<<";";
    }
    file<<"];"<<endl;
    
    file.close();
}

void output_statistics(string stat_file,real precision, real mean_int, int f_evals, int g_evals){
		// Quick fix, do not output sol( files grow too big)

    ofstream file;
    //if(nsol==1)
        //file.open(solfile.c_str(),std::ofstream::out);
    //else
    file.open(stat_file.c_str(),std::ofstream::out|std::ofstream::app);
    file.precision(16);
    
    // writing solution in matlab syntax
    file<<precision << ", " << mean_int << ", " << f_evals << ", " << g_evals << ";" << std::endl;
    
    file.close();
}
void output_mesh(string solfile)
{
	// Quick fix
	// don't output mesh
		return;
    ofstream file(solfile.c_str());
    file.precision(16);
    
    // writing mesh in matlab syntax
    file<<"x=[";
    for(int j=0;j<neqn;j++)
    {
        file<<((real)j)/(neqn-1.);
        if(j<neqn-1) file<<";";
    }
    file<<"];"<<endl;
    
    file.close();
}

void output_means(real mean, string mtfile)
{
    ofstream file(mtfile.c_str(),std::ofstream::out|std::ofstream::app);
    file.precision(16);
    
    file<<mean<<endl;
    
    file.close();
}
void output_evals_g(int evals, string evals_g_file){
    ofstream file(evals_g_file.c_str(),std::ofstream::out|std::ofstream::app);
    file.precision(16);
    
    file<<evals<<endl;
    
    file.close();
}

void output_evals_f(int evals, string evals_f_file){
    ofstream file(evals_f_file.c_str(),std::ofstream::out|std::ofstream::app);
    file.precision(16);
    
    file<<evals<<endl;
    
    file.close();
}

void square(real x[], real y[])
{
    for(int i=0;i<neqn;i++)
        y[i]=x[i]*x[i];
}

void exp(real x[], real y[])
{
    for(int i=0;i<neqn;i++)
        y[i]=exp(x[i]);
}

void sin(real x[], real y[])
{
    for(int i=0;i<neqn;i++)
        y[i]=sin(x[i]);
}

void arcsinh_2(real x[], real y[]){
	real tmp;
	for ( int i=0;i<neqn;i++){
		tmp=asinh(x[i]);
		y[i]=tmp*tmp;
	}
}

