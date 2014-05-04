/* 
 * File:   Ising5_2.h
 *
 */
//////////////////////////////////////////////////////////////////////////////
// Ising CLASS 
//////////////////////////////////////////////////////////////////////////////

#ifndef ISING5_2_H
#define	ISING5_2_H

#include <iostream>
#include <fstream>

#include <string>
#include <sstream>
#include <vector>
#include <cstdio>
#include <cmath>
#include <string>
#include <sstream>		// stringstream class
#include <cstring>		// Includes c_strings for file names
#include <cmath>
#include <iomanip>
#include "RNG.h"


using namespace std;

struct cluster_stat {
    vector<float> CF; //CF[x]=nb of occurences of the spin 
    //cluster of length x in a given state of the Ising chain 
    int cmax; //maximal length of the cluster 
} ;
// enum field used to choose kind of error calculated by ERROR function
enum ISING_ERROR_TYPE{
    ERROR_CHI,// calculate error of magnetic susceptibility 
    ERROR_CC  // calculate error of specific heat
};

    
class Ising {
public:
    Ising(){};
    Ising(int NN,double TT, double hh, double p, int maxGd, int maxTsep); //Constructor of the class
    virtual ~Ising();
    
    //***Ising spins chain and environmental parameters
    int N;              //chain length NN
    double T;           //temperature TT
    double h;           //magnetic field hh
    //double p;         //used by constructor to create initial state of the system
                        //p=0.5 correspond to high temperature limit, <S>=0
                        //p=1.0 correspond to low temperature limit, <S>=1
    
    vector<short> S;    //vector of Ising spins
    int maxGdist;       //limit the distance for which correlation functions Gs and Ge 
                        //will be calculated.    
    
//*** MEASURED OBSERVABLES AND MEASURING MEMBER FUNCTIONS ***
    //Measures are done for a given (instantaneous) state of the system
    //generated by Monte Carlo cycle. The instantaneous measures are 
    //collected and the mean value is taken after the last MC cycle.
    
    //ORDER PARAMETER
    double m;           //measured order parameter m=<S> (mean spin, magnetization)
    double order_parameter();   //measure m for a given spin configuration (instantaneous value)
    double mean_order_parameter();   //macroscopic observable
    

    //CORRELATION FUNCTIONS
    vector<double> S_correlation(int maxGdist,double OP); //function measuring instant spin-spin correlation
    vector<double> Gs;   //measured macroscopic correlation function Gs(spin separation)  
    vector<double> mean_S_corrrelation(string ); //macroscopic observable, returns Gs
    vector<double> Ge;   //measured correlation function Ge(local energy separation)  
    vector<double> E_correlation(int maxGdist,double E);  //function measuring instant energy-energy correlation
    vector<double> mean_E_correlation(int maxGdist,double E);  //function measuring macroscopic energy-energy correlation
    
    
    //ZERO FIELD SUSCEPTIBILITY
    double Chi(string totalFname);         //magnetic susceptibility calculated using order parameter fluctuation method
    double Chi( vector<vector<short> > SS );
    
    //SPECIFIC HEAT
    double CC();          //specific heat measured using energy fluctuation method calclutated from simulation
    double CC(string totalFname); // specific heat calculated from raw data file
    double CC(vector<vector<short> > SS); // specific heat calculated from raw data stored in memory
    
    //ENERGY
    virtual double E(); //measure energy per spin of a given spin configuration (lattice mean or instantaneous value) 
    double aveE(string totalFname); // energy from 
    
    //CLUSTER (DOMAIN) STATISTICS
    cluster_stat CFD;   //distribution CF of cluster length 
    cluster_stat cluster_freq(); //measure CF
    
    //TIME CORRELATION OF MONTE CARLO CYCLES
    double tau; //characteristic time of exponential decay of time autocorrelation function of magnetization
    vector<double> Gt; //autocorrelation function of magnetization
    int maxTimeseparation;
    vector<double> t_correlation(int maxTime, string ); //function measuring Gt
    
    //ERROR ANALYSIS, BOOTSTRAP METHOD
    double sigma; //standard deviation
    double ERROR(string totalFname,ISING_ERROR_TYPE error_type); //function returning sigma
    
//*** MONTE CARLO ENGINS, MC WRAPPER, SOME MC PARAMETERS 
    virtual vector<short> Metropolis_cycle();
    //void Wolff_cycle(int k);
    virtual void MC_simulation(int therm_t, int prod_t, int measure_f);
    int measure_f;      //frequency of measures, indicate number of MC cycles separating two consecutive measures
    int therm_t;        //thermalization time = number of initial cycles to thermalize (equilibrate) the initial state
    
    
    //Selected seed values for random number generators -- must belong to [1, 2^31 - 1]
    static long int seed1;// seed for random long integer generator
    static long int seed2;// seed for shuffling algorithm

    int nb_rejected;
    double rejected;
    

};

#endif	/* ISING5_2_H */

