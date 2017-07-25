#include <vector>
#include <iostream>
#include <complex>
#include <gsl/gsl_sf_fermi_dirac.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include "parameters.h"
#include "matrix.h"

using namespace LibConstants;
using namespace std;
// Some math functions and distributions

namespace MatLib {
	
 	double norm2(vector<double> vecnorm);
	
    // Linspace function generates vector
	vector<double> linspace(double start, double end, int num_in){
  		vector<double> linspaced;
  		double num = static_cast<double>(num_in);
    	double delta = (end - start) / (num - 1);
    
  		if (num == 0) { return linspaced; }
  			if (num == 1){
      			linspaced.push_back(start);
      		return linspaced;
    	}
  	    	for(int i=0; i < num-1; ++i){
            	linspaced.push_back(start + delta * i);
       		}
       
        linspaced.push_back(end); // start and end
                                  // are exactly the same as the input
        return linspaced;
    }
 		
	// Fermi distribution
 	double FermiDist(double E_value, double mu_value){
 		     double FDist = 1/(1 + exp((E_value - mu_value)/kT));
 		//   double FDist = n0*log(1+exp((mu_value - E_value)/kT));
 		return FDist;
	}
	

    // Fermi Dirac order +1/2
 	double FermiDiracHalf(double value){
 		double FDir = gsl_sf_fermi_dirac_half(value);
 		return FDir;
	 }
	
 	// Fermi Dirac order -1/2
 	double FermiDiracMHalf(double value){
 		double FDir = gsl_sf_fermi_dirac_mhalf(value);
 		return FDir;
	 }
	 
	// ACOS inverse trigonometric 
	complex<double> ARCCOS(double value){
		complex<double> result;
		gsl_complex sol = gsl_complex_arccos_real(value);
		result.real()=GSL_REAL(sol);
		result.imag()=GSL_IMAG(sol);
	
		return result;
	}
	
	// Solution k of the dispersion relation
	complex<double> k_disp(double dx, double E_level, double U_side, double E_k){
		complex<double> disp = ARCCOS(1 - m_e*dx*dx*(E_level - E_C - E_k - U_side)/(h*hh))/dx;
		return disp;
	}
	
	// Quasi-Fermi level
	vector<double>  QuasiFermiLevel(vector<double> E_Density, vector<double> Potential) {
		        unsigned Vec_size=E_Density.size();
                vector<double> QuasiFermi(Vec_size,0); 
				vector<double> lhs(Vec_size,0), root(Vec_size,0);
				double val;

                for(unsigned i=0;i<Vec_size;i++){
                    lhs[i] = E_Density[i]/N_C;
                    val = pow(3*sqrt(pi)*lhs[i]/4, (double)2/3); 
					root[i] = log(lhs[i]) / (1 - pow(lhs[i], (double)2)) + val / 
					          (  1 + pow(0.24 + 1.08*val, (double)(-2)) );
				    					
                    QuasiFermi[i] = kT*root[i] + Potential[i];    
                }
		 return QuasiFermi; 
    }
	
	// Poisson solver
	double Poisson(vector<double>& E_Density, vector<double>& Doping, 
							vector<double>& Potential, vector<double>& QuasiFermi){
								
		unsigned dim = Potential.size();
		double z;
		Matrix<double> D2(dim,dim,-2,1);
		
		
		vector<double> D(dim,0), D2U(dim,0), dN(dim,0), dU(dim,0);
			// Boundary conditions (Neumann) 
				D2(0,1)=2; D2(dim-1,dim-2)=2;
			//D2*U vector
			    D2U[0]     = -2*Potential[0] + 2*Potential[1];
				D2U[dim-1] =  2*Potential[dim-2] - 2*Potential[dim-1];
			//////////////////////////////////////////////////////
			for (unsigned i=0; i<dim; i++){
				z      =  (QuasiFermi[i] - Potential[i])/kT; 
				D[i]   =  q*dx*dx/(eps*kT)*N_C*FermiDiracMHalf(z);
					if ( (i!=0) && (i!=dim-1) ){
						D2U[i] =  Potential[i-1] - 2*Potential[i] + Potential[i+1];
					}
				dN[i]    = -(q*dx*dx/eps)*(E_Density[i]-Doping[i])-D2U[i]; 
				D2(i,i) -= D[i];
			}
			
			              dU = D2/dN;
						  
		// Add correction
			for (unsigned i=0; i<dim; i++){
				Potential[i] += dU[i];
			}
			
			z = norm2(dU)/norm2(Potential); 
			
		return z;
	}
	
	
	
	// Display vestor
	void show(vector<double> vector_input){
		for(unsigned i=0; i<vector_input.size(); i++){
			cout << vector_input[i] << " ";
		} cout << endl;
	}
	
	// l2-norm
	double norm2(vector<double> vecnorm){
		double norm = 0.0;

         for (int i = 0; i < vecnorm.size(); ++i) {
                 norm += vecnorm[i] * vecnorm[i];
         }

       return sqrt(norm);
	}
	

};
