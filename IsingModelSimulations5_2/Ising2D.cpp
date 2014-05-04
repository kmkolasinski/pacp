#include "Ising2D.h"


//////////////////////////////////////////////////////////////////////////////
//  section A             CONSTRUCTOR OF OBJECT "ISING"                     //          
//////////////////////////////////////////////////////////////////////////////


Ising2D::Ising2D(){

}

Ising2D::Ising2D(int NN, double TT, double hh, double p, int maxGd, int maxTsep) {
    //CONSTRUCTOR FOR 1 DIMENSIONAL ISING CHAIN
    N = NN*NN; //set length of the square box
    T = TT;    //set temperature
    h = hh;    //set magnetic field
    maxGdist = maxGd;
    maxTimeseparation = maxTsep;

    //INITIALISATIONS 
    //***ADD OTHER VARIABLES AND INITIALISATIONS WHEN NECESSARY 

    m = 0; //set to zero initial value of magnetization

    //Creation of initial state of the chain
    RNG::initializeRNG(seed1); // initializes integer random number generator
    for (int i = 0; i < N; i++) {
        if (RNG::drandom(seed2) <= p)
            S.push_back(1);
        else
            S.push_back(-1);
    }

    //Set to zero initial value of spin-spin and energy-energy correlation function
    for (int i = 0; i < maxGdist; i++) {
        Gs.push_back(0);
        Ge.push_back(0);
    }

    //cout<<"  S "<<S[0];


    //Set to zero initial value of time auto-correlation function of magnetization
    for (int i = 0; i < maxTimeseparation; i++) {
        Gt.push_back(0);
    }

    //Set to zero initial value of cluster statistic vector
    for (int i = 0; i < N; i++) {
        CFD.CF.push_back(0);
        CFD.cmax = 1;
    }

    
    
    
}

Ising2D::~Ising2D() {
}





//////////////////////////////////////////////////////////////////////////////
// section B                 MONTE CARLO                                    //
//////////////////////////////////////////////////////////////////////////////

vector<short> Ising2D::Metropolis_cycle() {
    /* Domain : Monte Carlo simulation in 2D
     * Perform Metropolis update of N randomly chosen lattice sites 
     * i.e. one MC cycle  (one sweep of the lattice)*/

    int site, siteL, siteR , siteU, siteD; //site=randomly chosen site, added up and down
    //siteL its left interacting neighbor,
    //siteR its right interacting neighbor.
    double deltaE;

    for (int i = 0; i < N; i++) {
        site = (int) (N * RNG::drandom(seed2));
        siteL = (site - 1 + N) % N; //periodic boundary conditions
        siteR = (site + 1) % N; //periodic boundary conditions
        deltaE = 2 * S[site]*(S[siteL] + S[siteR]); //change of energy if spin
        //at site "site" is flipped (units : beta*J -> 1/T, J -> 1, where J is coupling constant)    
        if (deltaE <= 0) {
            S[site] = -S[site];
        } else {
            if (RNG::drandom(seed2) < exp(-deltaE / T)) {
                S[site] = -S[site];
            } else
                nb_rejected++;
        }
    }
    return S; //The state of the system after one lattice sweep is returned
}