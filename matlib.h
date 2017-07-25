#ifndef _MAT_LIB_H
#define _MAT_LIB_H

using namespace std;
// Some math functions and distributions
namespace MatLib {
 	// Linspace function generates vector
	vector<double> linspace(double start, double end, int num_in);
 	// Fermi distribution
 	double FermiDist(double E_value, double mu_value);
	// Fermi Dirac order +1/2
 	double FermiDiracHalf(double value);
	// Fermi Dirac order -1/2
 	double FermiDiracMHalf(double value);
 	// ACOS inverse trigonometric 
	complex<double> ARCCOS(double value);
	// Solution k of the dispersion relation
	complex<double> k_disp(double dx, double E_level, double U_side, double E_k);
	// Quasi-Fermi level
	vector<double>  QuasiFermiLevel(vector<double> E_Density, vector<double> Potential);
	// Poisson solver
	double Poisson(vector<double>& E_Density, vector<double>& Doping,
							vector<double>& Potential, vector<double>& QuasiFermi);
	// Norm L2 vector
	double norm2(vector<double> vecnorm);
	//Show vector
	void show(vector<double> vector_input);
};
#endif