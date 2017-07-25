#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <complex>
#include <essl.h>

#include <mpi.h>

#include "matrix.h"
#include "parameters.h"
#include "matlib.h"



#define _i cmplx(0,1)
#define PoissonConvergence 1e-5


typedef complex<double> cmplx;
// Namespaces 
using namespace LibConstants;
using namespace MatLib;

// Dyson Equation
	void DysonEquation(Matrix<cmplx>& G_R, double k_point, double E, const Matrix<cmplx>& Hamiltonian, 
						       vector<double> Potential);
	
// Keldysh Equation
	void KeldyshEquation(Matrix<cmplx>& G_les, Matrix<cmplx>& BSE_les, double k_point, double E, 
							Matrix<cmplx>& G_Rin, vector<double> Potential, double Fermi_left, double Fermi_right);

int main(int argc, char **argv) {
	
	// Energy space
	vector<double> Energy=linspace(Emin,Emax,NE);
	double DE = Energy[1] - Energy[0];
	
	// Momentum space
	vector<double> k_t=linspace(0,Kmax,NK);
	double DK = k_t[1] - k_t[0];
	
	
	// Fermi level and bias Vp [Volts]
	double Vp = 0.2;
	double Fermi_L = mu, Fermi_R = mu-Vp;
	
	// Define potential, electron density, doping concentration
	const double  nd = N_C*FermiDiracHalf((mu-E_C)/kT);    // reference doping level (outer regions)
	vector<double> N_D(Nc,0.1*nd); 
	// Doping profile (add sides)
	N_D.insert(N_D.begin(),Ns,nd);
    N_D.insert(N_D.end(),Ns,nd);
	
	//Poisson epsilon
    double PoissonEps;
    
	
	// NEGF Matrices initialized
	const Matrix<cmplx> H(Nx, Nx, E_C + 2*t, -t), Zero(Nx,Nx,0);
	
	Matrix<cmplx> G_R(Nx,Nx,0), G_les(Nx,Nx,0);
	Matrix<cmplx> BSE_les(Nx,Nx,0);
	Matrix<cmplx> G_RT(Nx,Nx,0), G_lesT(Nx,Nx,0);
	
	// Observables
	vector<double> 	U(Nx,0), 
					QuasiFermi(Nx,0),	
					ElectronDensityTotal_E(Nx,0), 	   ElectronDensityTotal(Nx,0);
					
	double 			ElectronDensitySpectral_E[Nx],     ElectronDensitySpectral[NE][Nx], 
					CurrentDensitySpectral_E[NE][Nx],  CurrentDensitySpectral[NE][Nx],
					LDOS_E[Nx],                        LDOS[NE][Nx]; 
	//////////////////////////////////////////////////////////////////////////////////
	
	
	
	
	//MPI time control
	double starttime, endtime; 	
    
    
	MPI_Init(NULL,NULL);
	int pid;
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	int mpi_size;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
		
	////////////////////////////
	
	
	
	int Cycle=0;
	
	starttime = MPI_Wtime();	
	
	do{	
	
	
	if (pid > 0){
		
		G_lesT = Zero;
		G_RT   = Zero;
		
		for(unsigned k=0; k<NK; ++k){
					
			DysonEquation(G_R, k_t[k], Energy[pid-1], H, U);
			
			KeldyshEquation(G_les, BSE_les, k_t[k], Energy[pid-1], G_R, U, Fermi_L, Fermi_R);
			
			G_lesT   += G_les*(DK*k_t[k]/(2*pi));
			
			G_RT     += G_R*(DK*k_t[k]/(-2*pi*dx));
			
		}
		
		
		// Density contribution from one energy nE
		for (int j=0; j<Nx; j++){
		       ElectronDensityTotal_E[j]  	= G_lesT(j,j).imag()*2*DE/(2*pi*dx);
			   ElectronDensitySpectral_E[j] = G_lesT(j,j).imag()*2/dx;
			   LDOS_E[j] 					= G_RT(j,j).imag()*2/(pi*dx);
		}
		
		
		    MPI_Send(&ElectronDensitySpectral_E, Nx, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);  
			MPI_Send(&LDOS_E, Nx, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		
	}
	
		
		MPI_Reduce(&ElectronDensityTotal_E[0], &ElectronDensityTotal[0], 
						ElectronDensityTotal.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
	
		if (pid == 0) {
			Cycle++;
			
			for (unsigned np=0; np<NE; np++){
				MPI_Recv(&ElectronDensitySpectral[np], Nx, MPI_DOUBLE, np+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&LDOS[np], Nx, MPI_DOUBLE, np+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}			  
			
			
			// Quasi-Fermi level
			QuasiFermi=QuasiFermiLevel(ElectronDensityTotal, U);
					
			
		    // Poisson solution 
			PoissonEps = Poisson(ElectronDensityTotal, N_D, U, QuasiFermi);
			
					
			cout << "Cycle# " <<  Cycle  << " Convergence " << PoissonEps << endl;
			
			
			
		
			
			if (PoissonEps<PoissonConvergence){
				
				ofstream Potential_out, EDT_out, EDS_out, LDOS_out;
				Potential_out.open 	("Potential.dat");
				EDT_out.open 		("ElectronDensityTotal.dat");
				EDS_out.open 		("ElectronDensitySpectral.dat");
				LDOS_out.open 		("LDOS.dat");
				
				endtime = MPI_Wtime();
					
					cout <<"=========================================" << endl;;
                    cout <<"Wall-clock time: " << (endtime-starttime)  << " seconds" << endl;
                    cout <<"=========================================" << endl;;
					
					for (int j=0; j<Nx; j++){
						Potential_out << j << " " << U[j] << endl;
						EDT_out << j << " " << ElectronDensityTotal[j] << " " << N_D[j] << endl;
					}
					
					for (int i=0; i<NE; i++){
						for (int j=0; j<Nx ;j++){
								EDS_out 	<< 	ElectronDensitySpectral[i][j]	 << " ";
								LDOS_out 	<<	LDOS[i][j]						 <<	" ";
						}
						        EDS_out 	<<	endl;
								LDOS_out 	<<	endl;
					}
					
					Potential_out.close();
					EDT_out.close();
					EDS_out.close();
					LDOS_out.close();
					
					
			}
			
			 
			 
			
        }
		
		
		MPI_Bcast(&U[0], U.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&PoissonEps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		

			
	}while(PoissonEps>PoissonConvergence);	
		
	MPI_Finalize();  
	
	
			
  return 0;
}





// Dyson Equation
	
	void DysonEquation(Matrix<cmplx>& G_R, double k_point, double E, const Matrix<cmplx>& Hamiltonian, 
						       vector<double> Potential)
	{
		
		cmplx BSE_L, BSE_R;
		double E_K = h*hh*k_point*k_point/(2*m_e);
		
		// Initialize Green's function
		G_R.init_zero();
		
		// Boundary self-energies		
		BSE_L=-t*exp(_i*k_disp(dx, E , Potential[0], E_K)*dx);
   		BSE_R=-t*exp(_i*k_disp(dx, E, Potential[Nx-1], E_K)*dx);
		
		// Dyson equation
		G_R -= Hamiltonian; G_R(0,0) -= BSE_L; G_R(Nx-1,Nx-1) -= BSE_R;  
		
		for (int j=0; j<Nx; j++){ 
		      G_R(j,j) += E - E_K - Potential[j];
		}
		
		G_R.inverse();	
		// End Dyson equation
				
	} 
	
// Keldysh Equation
	void KeldyshEquation(Matrix<cmplx>& G_les, Matrix<cmplx>& BSE_les, double k_point, double E, 
							Matrix<cmplx>& G_Rin, vector<double> Potential, double Fermi_left, double Fermi_right)
	{
		
		// Transverse momentum
		double E_K = h*hh*k_point*k_point/(2*m_e); 
		// Boundary self-energier (only corner elements)
		cmplx BSE_left, BSE_right;
		
		BSE_left =-t*exp(_i*k_disp(dx, E , Potential[0], E_K)*dx);
   		BSE_right=-t*exp(_i*k_disp(dx, E, Potential[Nx-1], E_K)*dx);
		// Boundary self energy corners prescribed
		BSE_les(0,0)       = -2.0*_i*imag(BSE_left)*FermiDist(E, Fermi_left);
		BSE_les(Nx-1,Nx-1) = -2.0*_i*imag(BSE_right)*FermiDist(E, Fermi_right);
		// Equation itself
		G_les = G_Rin; G_les *= BSE_les; G_les *= G_Rin.dagger();
		// Return same G_R as required
		G_Rin.dagger();

	} 

