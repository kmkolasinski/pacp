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
    nx = NN;
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

    int site,site_i,site_j, siteL, siteR , siteU, siteD; //site=randomly chosen site, added up and down
    //siteL its left interacting neighbor,
    //siteR its right interacting neighbor.
    //siteU its up interacting neighbor.
    //siteD its down interacting neighbor.
    double deltaE;    
    for (int i = 0; i < N; i++) {
        site_i = (int) (nx * RNG::drandom(seed2));
        site_j = (int) (nx * RNG::drandom(seed2));
        site  = To2D(site_i,site_j);
        
        // periodic representation
        siteL = To2D((site_i - 1 + nx)%nx,site_j); 
        siteR = To2D((site_i + 1     )%nx,site_j); 
        siteD = To2D(site_i,(site_j- 1 + nx)%nx); 
        siteU = To2D(site_i,(site_j+1      )%nx); 
        
        // helicoidal 
        //siteL = (To2D(site_i - 1,site_j)  + N) % N; 
        //siteR = (To2D(site_i + 1,site_j)     ) % N; 
        //siteD = (To2D(site_i ,site_j-1)   + N) % N; 
        //siteU = (To2D(site_i ,site_j+1)      ) % N; 
             
        
        deltaE = 2 * S[site]*(S[siteU]+S[siteD]+S[siteL] + S[siteR]); //change of energy if spin
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


void Ising2D::MC_simulation(int therm_t, int prod_t, int measure_f) {
    /* Domain : Metropolis simulation of Ising 2d model
     * Perform desired number of MC cycles
     * therm_t = thermalization cycles
     * prod_t  = production cycles
     * Store raw results (states of the chain) in the disk file. */

    RNG::initializeRNG(seed1); // initializes integer random number generator
    
    cout << "----------------------------------------" << endl;
    cout << "MC simulation start..." << endl;
    cout << "Thermalization time:" <<  therm_t    << endl;
    cout << "Production     time:" <<  prod_t     << endl;
    cout << "Measure frequency  :" <<  measure_f  << endl;
    
    
    //*** INITIAL MC CYCLES ***
    //Thermalization of initial non-equilibrium state.
    for (int t = 0; t < therm_t; t++) {
        Metropolis_cycle();
    }
    cout << "Thermalization done..." << endl;
    //SOME LOCAL VARIABLES; ADD NEW WHEN NECESSARY
    nb_rejected = 0; //Provides information about the efficiency of sumulation

    //STORING MC DATA IN A FILE
    string path, fname, totalFname;
    //path: use the path to the project folder on your laptop
    path = "";
    fname = "MCdata";
    totalFname = path + fname;
    fstream DATA(totalFname.c_str(), ios::out); //file is open for writing
    DATA.setf(ios::showpos);

    //*** MONTE CARLO PRODUCTION CYCLES ***
    //Turns on prod_t times MC engine and accumulate results for observables.
    vector<short> s(N, 0);
    int progress = 0;
    for (int t = 0; t < prod_t; t++) {
        if (!(t % measure_f)) {
            s = Metropolis_cycle();
            for (int i = 0; i < N; i++)
                DATA << s[i] << " ";
        } else
            Metropolis_cycle();
        /*
         Please make copy of this file (Ising5_2.cpp)
         */
        DATA << endl;
        if(t%(prod_t/10) == 0){
            cout << progress << "%\t";
            cout.flush();
            progress += 10;
        }
    }
    cout << endl;
    DATA.close();
}


//////////////////////////////////////////////////////////////////////////////
// section C                  DATA ANALYSIS                                 //          
//////////////////////////////////////////////////////////////////////////////

double Ising2D::E(){    
    double E;
    E=0;
    for(int i=0;i<nx;i++){
    for(int j=0;j<nx;j++){        
        E+= -S[To2D(i,j)]*( S[To2D((i+1)%nx,j)] + S[To2D(i,(j+1)%nx)] ) + h*S[To2D(i,j)]; 
    }}
    return E/(N);
}