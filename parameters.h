#ifndef _LIB_CONSTANTS_H
#define _LIB_CONSTANTS_H


namespace LibConstants{

    // Fundamental Physical Constants
 const double  h     =   6.582119514e-16;   // Plank's constant [eV*s]
 const double  hh    =   1.0545718e-34;     // Plank's constant [J*s]
 const double  m0    =   9.10938356e-31;    // electron mass
 const double  q     =   1.60217662e-19;    // elementary charge [C]
 const double  kT    =   0.0256;            // Boltzmann constant x temperature T=300K [eV]
 const double  pi    =   3.14159265358979;  // Pi

    // GaAs Physical Constants
 const double m_e    = 0.067*9.10938356e-31;    // effective mass [kg]
 const double eps    = 10.5*8.85e-12;           // permittivity [F/m]
 const double E_C    = 1.420;                   // band edge energy [eV]
 
   // Structure / discretization
 const  unsigned  Ns   = 15;                     // Number of points in doped section
 const  unsigned  Nc   = 90;                     // Number of points in central section
 const  unsigned  Nx   = Ns + Nc + Ns;           // Total number of points 
 const  double    dx   = 1e-9;                   // Step in space
   // Energy discretization
 const  unsigned  NE   = 1000;                   // Number of energy points
 const  double    Emin = E_C - 0.22;             // Energy range min
 const  double    Emax = E_C + 0.58;             // Energy range max
  // Fermi level
 const  double    mu   = E_C+0.4;                // Fermi level
  // Momentum discretization
 const  unsigned  NK    = 100;                     // Transverse momentum points
 const  double    Kmax  = 1.5/1e-9;               // Maximum k-value [nm]^(-1)
 
   
   // Other physical parameters (derived)
 const double  n0    = 2*m_e*kT*q/(2*pi*hh*hh);      
 const double  t     = h*hh/(2*m_e*dx*dx);         // Hopping parameter for Effective mass Hamiltonian
 const double  N_C   = 2*std::pow(0.5*n0,1.5);     // Effective density of states in the conduction band
 
 
 
}

#endif
