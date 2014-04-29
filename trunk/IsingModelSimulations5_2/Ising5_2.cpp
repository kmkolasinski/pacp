/* 
 * File:   Ising5_2.cpp
 */
//////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION OF Ising CLASS MEMBER FUNCTION
//This file contains four sections :
// section A    CONSTRUCTOR OF OBJECT "ISING"                                         
// section B    MONTE CARLO                                                 
// section C    DATA ANALYSIS                                                         
// section D    DATA ANALYSIS: EXAMPLE. 
// Sections A, B and D should stay as they are
// Section B contain Metropolis engine (Wolff engine will come soon) 
// The functions of section B will be used in "main()" function to produce
// time series of states of the chain created by Metropolis algorithm.
// The series will be stored in disk file and used by functions of section C.
// The functions of section C has to be created by teams as indicated, however,
// hungry people can do more! No restrictions for this. Good ideas are welcome.
// Just take into account that project has to be developed according to standards 
// of collaborative work. We are its contributors. 
//////////////////////////////////////////////////////////////////////////////

/*** PLEASE MAKE COPY OF THIS FILE. USE IT ON YOUR LOCAL COMPUTER 
 *** TO WRITE, ADD,  AND TEST FUNCTIONS.
 *** 
 *** HANG BACK THE FUNCTIONS YOU HAVE CREATED ON OUR EXCHANGE BLACKBOARD.*/

#include <iostream>
#include <fstream>
#include "stdio.h"
#include "math.h"
#include <string>
#include <sstream>		// stringstream class
#include <cstring>		// Includes c_strings for file names
#include <vector>
#include <cmath>
#include <iomanip>
#include "RNG.h"
//#include<stack> //LIFO (Last In First Out) container for Wolff MC algorithm

#include "Ising5_2.h"
//#include "Tools.h"



using namespace std;

//Selected seed values for random number generators -- must belong to [1, 2^31 - 1]
long int seed1 = 46723; // seed for random long integer generator
long int seed2 = 2593459; // seed for shuffling algorithm

vector<vector<short> > readDATAtoVectors(string totalFname, int N); //forward declaration (Tools)




//////////////////////////////////////////////////////////////////////////////
//  section A             CONSTRUCTOR OF OBJECT "ISING"                     //          
//////////////////////////////////////////////////////////////////////////////

Ising::Ising(int NN, double TT, double hh, double p, int maxGd, int maxTsep) {
    //CONSTRUCTOR FOR 1 DIMENSIONAL ISING CHAIN
    N = NN; //set length of the chain
    T = TT; //set temperature
    h = hh; //set magnetic field
    maxGdist = maxGd;
    maxTimeseparation = maxTsep;

    //INITIALISATIONS 
    //***ADD OTHER VARIABLES AND INITIALISATIONS WHEN NECESSARY 

    m = 0; //set to zero initial value of magnetization

    //Creation of initial state of the chain
    initializeRNG(seed1); // initializes integer random number generator
    for (int i = 0; i < N; i++) {
        if (drandom(seed2) <= p)
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

Ising::~Ising() {
}






//////////////////////////////////////////////////////////////////////////////
// section B                 MONTE CARLO                                    //
//////////////////////////////////////////////////////////////////////////////

vector<short> Ising::Metropolis_cycle() {
    /* Domain : Monte Carlo simulation
     * Perform Metropolis update of N randomly chosen lattice sites 
     * i.e. one MC cycle  (one sweep of the lattice)*/

    int site, siteL, siteR; //site=randomly chosen site, 
    //siteL its left interacting neighbor,
    //siteR its right interacting neighbor.
    double deltaE;

    for (int i = 0; i < N; i++) {
        site = (int) (N * drandom(seed2));
        siteL = (site - 1 + N) % N; //periodic boundary conditions
        siteR = (site + 1) % N; //periodic boundary conditions
        deltaE = 2 * S[site]*(S[siteL] + S[siteR]); //change of energy if spin
        //at site "site" is flipped (units : beta*J -> 1/T, J -> 1, where J is coupling constant)    
        if (deltaE <= 0) {
            S[site] = -S[site];
        } else {
            if (drandom(seed2) < exp(-deltaE / T)) {
                S[site] = -S[site];
            } else
                nb_rejected++;
        }
    }
    return S; //The state of the system after one lattice sweep is returned
}

void Ising::MC_simulation(int therm_t, int prod_t, int measure_f) {
    /* Domain : Metropolis simulation of Ising 1d model
     * Perform desired number of MC cycles
     * therm_t = thermalization cycles
     * prod_t  = production cycles
     * Store raw results (states of the chain) in the disk file. */

    initializeRNG(seed1); // initializes integer random number generator

    //*** INITIAL MC CYCLES ***
    //Thermalization of initial non-equilibrium state.
    for (int t = 0; t < therm_t; t++) {
        Metropolis_cycle();
    }

    //SOME LOCAL VARIABLES; ADD NEW WHEN NECESSARY
    nb_rejected = 0; //Provides information about the efficiency of sumulation

    //STORING MC DATA IN A FILE
    string path, fname, totalFname;
    //path: use the path to the project folder on your laptop
    path = "/home/tomek/Desktop/NetBeansProjects/IsingModelSimulations5_2/";
    fname = "MCdata";
    totalFname = path + fname;
    fstream DATA(totalFname.c_str(), ios::out); //file is open for writing
    DATA.setf(ios::showpos);

    //*** MONTE CARLO PRODUCTION CYCLES ***
    //Turns on prod_t times MC engine and accumulate results for observables.
    vector<short> s(N, 0);
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
    }
    DATA.close();
}


//void Ising::Wolff_cycle(int k) {
//    /* Domain : Monte Carlo simulation
//     * Perform Wolff updates of cluster of spins 
//     */
//     /*
//      * To be completed later
//      */
//}







//////////////////////////////////////////////////////////////////////////////
// section C                  DATA ANALYSIS                                 //          
//////////////////////////////////////////////////////////////////////////////

double Ising::E(){
    double E;
    E=0;
    for(int i=0;i<N;i++){
        if(i==(N-1)){
            E+=-S[i]*S[0];
        }else{
            E+=-S[i]*S[i+1];
        }
    }
    return E/(N);
}

vector<double> Ising::E_correlation(int maxGdist, double E) {
    /*Domain: data analysis, calculation of physical quantities
     * Calculate instantaneous connected correlation function
     * E represents instantenous (i.e. lattice) mean value of energy. 
     *         Ge(x)= <E(0)E(x)> - <E>^2
     *** PLEASE COMPLETE THIS COMMENT
     */

    double E2 = E*E;
    vector<double> corr;
    /*
     * Panowie Kaczmarczyk i Jedraczka we wspolpracy z 
     * pania Chadrian z kolezanka.
     * 
     */
    return corr;

}


vector<double> mean_E_correlation(int maxGdist,double E){
    /*
     * Panowie Kaczmarczyk i Jedraczka we wspolpracy z 
     * pania Chadrian z kolezanka.
     * 
     */
}


double Ising::CC() {
    /*
     * returns value of the specific heat according to relation:
     * C=bheta/T <E^2>-<E>^2
     * Function is using Metropolis_cycle() to sweep through latice between 
     * measurements. The value if acquired with frequency declared as measure_f
     * variable in main function.
     * this function needs E() to work properly.
     */
    double meanE,meanEsq;
    int n=0;
    for (int t = 0; t < 1000; t++) {
        if (!(t % 5)){
            meanE+=E();
            meanEsq+=meanE*meanE;
            Metropolis_cycle();
            n++;
        }else{
            Metropolis_cycle();
        }
    }
    meanE/=n;
    meanE*=meanE;
    meanEsq/=n;
    return (meanEsq-meanE)/T;
}



double Ising::Chi(string totalFname) {
    // Calculates the magnetic susceptibility of the system
    // which is stored in file totalFname
    //parameter: totalFname - name of file
    
    //Raw simulation data flows from disk to vector container
    vector<vector<short> > SS = readDATAtoVectors(totalFname, N);

    return Chi(SS);
}

double Ising::Chi( vector<vector<short> > SS ) {
    // Calculates the magnetic susceptibility of the system
    // which is stored in vector SS
    //parameter: SS - Raw simulation data
    
    //number of time samples
    int t_steps=SS.size();
    //auxiliary variables
    double M_temp, aM, aM2;
    
    aM=0.; aM2=0.;
    for(int j=0;j<t_steps;j++){
        M_temp=0.;
        //calculation of magnetisation M
        for(int i=0;i<N;i++){
            M_temp+=SS[j][i];
        }
        //calculation of <M> and <M^2>
        aM+=M_temp;
        aM2+=M_temp*M_temp;
    }
    //normalization by number of measures and sites
    aM=aM/t_steps/N; 
    aM2=aM2/t_steps/N;
    
    return (aM2-aM*aM)/T; //Chi
}

vector<double> Ising::t_correlation(int maxTime, string totalFname) {
    /*Domain: data analysis, calculation of MC characteristic
     * PLEASE WRITE HERE COMMENTS */

    vector<double> mt(maxTime, 0.);
    //

    /* RECIPE FOR THE CALCULATION OF THE CORRELATION TIME
     * The autocorrelation Psi defined as Psi(t) = \int{dt'[m(t')m(t+t')]-<m>^2},
     * where m(t) is instantaneous value of magnetization and <m> is the mean value.
     * The integration runs from zero to the total time of simulation. 
     * The time displacement "t" takes values from 0 to maxTimeseparation,
     * where maxTimeseparation has to be established in the following way:
     * Psi values between maxTimeseparation/2 and maxTimeseparation should
     * stabilize and oscilate around zero, while for times t less then maxTimeseparation/2
     * exponential decay should be observed. The characteristic time tau of  
     * exponential decay has to be find by further graphical analysis of 
     * log(psi(t)) as function of t. Independently the check may be done by 
     * integration of psi(t). Psi(t)~exp(-t/tau) implies that integral of Psi(t)/Psi(0)
     * from 0 to infinity gives 1/tau.
     * 
     */

     /*
     * Pan Mateusz Dyndal 
     */
     
   //Raw simulation data flows from disk to vector container
    vector<vector<short> > SS = readDATAtoVectors(totalFname, N);
    
    //Declaration of iterators
    vector<vector<short> >::iterator t;
    vector<short>::iterator pos;
    
    int i=0;

    for (t = SS.begin(); t != SS.end(); t++) { //loop over time "t" (time unity in MC cycles = prod_t/measure_f)
        for (pos = t->begin(); pos != t->end(); pos++) { //loop over chain sites "pos" at fixed time
            mt[i] += (double)(*pos);
        }
        mt[i] /= (double)N;
        i++;
        //if(i==maxTimeseparation) break;
    }


    cout << "\n--------------------";

    for (int j = 0; j < maxTimeseparation; j++)
        cout << "\nmt[" << j << "]=" << mt[j];

    cout << "\n--------------------";



    double meanm =0.;
    for (int j = 0; j < maxTime; j++)
        meanm += mt[j];
    meanm/=(double)maxTime;


/////////////////// MAIN LOOP ///////////////////////////////////
    for (int j = 0; j < maxTimeseparation; j++){
    
      for (int k = 0; k < maxTime; k++){  
       
       int index = j + k;
       if ((j + k) >= maxTime) index -= maxTime;
       
       Gt[j] += (mt[k]*mt[index] - meanm*meanm);
       
      }
    
    }
/////////////////////////////////////////////////////////////////
    
    double Gt_int = 0.;
    fstream DATA("Gt.dat", ios::out);
    for (int j = 0; j < maxTimeseparation; j++){
        cout << "\nGt[" << j << "]=" << Gt[j];
        DATA << Gt[j] << "\n";
        Gt_int += Gt[j];
    }
    cout << "\n--------------------";

    cout<<"i = "<<i<<endl;
    cout<<"meanm = "<<meanm<<endl;
    cout<<"Gt integral = "<<Gt_int<<endl;
    Gt_int /=Gt[0];
    cout<<"tau = "<<Gt_int<<endl;
    Gt_int = 1./Gt_int;
    cout<<"1/tau = "<<Gt_int<<endl;
    
    return Gt;

}

//function returning sigma
double Ising::ERROR(string totalFname) {
   /* BOOTSTRAP RECIPE OF ERROR CALCULATION
    * Let n be a number of elements in dataSet.
    * 1. Pick at random n elements (with returns).
    * 2. Calculate an observable C using n elements
    *    created in 1.
    * 3. Repeat 1 and 2 m times. This gives a series
    *    C1, C2, ..., Cm of estimations of C.
    * 4. Use this series to calculate the standard
    *    deviation sigma in the following way:
    *    sigma= sqrt( <C^2> - <C>^2 )
    */
   
//Raw simulation data flows from disk to vector container
   vector<vector<short> > SSdata = readDATAtoVectors(totalFname, N);
   vector<vector<short> > SS = SSdata;
   int n = SS.size();
    
   double chi ;
   double sum2=0, sum1=0;
   int m=100;
   
   for(int j=0;j<m;j++){
      
        // we choose random set (with returns)
        for( int i = 0 ; i < n ; i++ ){
           int time_step = (int) ( n * drandom(seed2));
           SS[i] =  SSdata[time_step];              
        }
   
        chi = Chi(SS);

        sum2 += chi*chi;
        sum1 += chi;
   }
   sum2 /= m;
   sum1 /= m;
   
   double sigma = sqrt( sum2 - sum1*sum1 );
   
    return sigma;


}

cluster_stat Ising::cluster_freq() {
    /* Domain: data analysis. 
     * Goal: calculation of the cluster size distribution.
     * For 1-dimensional system T=0 corresponds to criticality
     * When approaching critical state the system develops large clusters of 
     * spins pointing in the same directions. Such clusters are named domains.
     * When going closer and closer to critical temperature (T=0 in 1d case)
     * clusters are not only growing in size but their distribution flatten
     * meaning that domains of any size are present close to criticality.
     * This function is designed to check this behavior. It calculates
     * the domain length instant distribution. Function returns structure
     * "cluster_stat" containing two members: 
     * - vector CF: CF[L] contain the frequency of domain of size L
     * - cmax : contain the size of the largest domain in the distribution
     */

    //***Please use local variables as defined below
    vector<float> tt(N, 0); //tt is created to initialize the vector component
    //of the cls instance of the cluster_stat structure below
    cluster_stat cls = {tt, 1}; //container for results
    int pos=0; //initial current position in the chain
    int sg = S[pos]; //initial sign of current position (i.e. spin orientation)
    int cl_L = 1; //initial length of the current cluster
    int max=0;//temporary variable used in the for loop
    int amax=0;//absolute maximum value of domain size
    float mmax=0;//mean maximum value of domain size
    int stat=1;//statistic of spin domain orientation in fixed temperature
    for(int i=0;i<stat;i++){
        pos=0;
        while(pos<N){
            if(pos<(N-1)){///periodic boundary conditions
                if(S[pos++]==sg)cl_L++;
                else{
                    tt[cl_L]+=1;
                    if(max<cl_L) max=cl_L;
                    cl_L=1;
                }
            }else{
                if(S[pos]==S[++pos])cl_L++;
                else{
                    tt[cl_L]+=1;
                    if(max<cl_L) max=cl_L;
                    cl_L=1;
                }
            }
        }
        if(amax<max) amax=max;
        mmax+=max;
        max=0;
        Metropolis_cycle();
    }
    for(int i=0;i<N;i++){
        tt[i]/=stat;
    }
    cls.CF=tt;
    cls.cmax=amax;
    mmax/=stat;

    //If the first and the last spin of the chain point in the same direction
    //they belong to the same domain (we consider chains with periodic boundary conditions)

    /*
     * 
     Panowie Kaczmarczyk i Jedraczka
     * 
     */

    return cls;

    //TESTS        
    //        Display(S);
    //        for (int i = 0; i <= cls.cmax; i++) {
    //           if (cls.CF[i]) cout<<"\n nb of blocks of length "<<i<<" = "<< cls.CF[i];
    //        }
    //        int sum=0;
    //        for (int i = 0; i <= cls.cmax; i++) {
    //           sum+=i*cls.CF[i];
    //        }
    //        cout<<"\n test of blocks statistic: sum = "<<sum;
}






//////////////////////////////////////////////////////////////////////////////
// section D                DATA ANALYSIS: EXAMPLE                          //          
//////////////////////////////////////////////////////////////////////////////

double Ising::order_parameter() {
    /*Domain: data analysis, calculation of physical quantities
      Calculation of instantaneous lattice mean value <S>
       of the spin, (i.e. order parameter).
     */
    double OP; //local name for order parameter
    for (int i = 0; i < N; i++) {
        OP += S[i];
    }
    return OP /= N;
}

vector<double> Ising::S_correlation(int maxGdist, double OP) {
    /* Domain : calculation of physical quantities.
     * Calculate instantaneous connected correlation function 
     *           G(x)= <S(0)S(x)>-<S(0)><S(x)>
     * for homogeneous system. Works on Ising data member S(i)
     * i=1..N, representing an instantaneous state of the 
     * chain. Input parameter OP (Order Parameter) represents
     *  i.e. instantaneous (lattice) mean spin value <S> calculated
     * in "order_parameter" function. Function "S_correlation" returns a vector 
     * containing instantaneous values of correlations ranging from 0
     * to maxGdist. For small systems maxGdist should be smaller then N/2.*/

    double OP2 = OP*OP;
    vector<double> corr(maxGdist, 0);
    for (int dist = 0; dist < maxGdist; dist++) {
        for (int pos = 0; pos < N; pos++) {
            corr[dist] += S[pos] * S[(pos + dist) % N];
        }
        corr[dist] /= N;
        corr[dist] -= OP2;
    }
    return corr;

}

// vector<vector<short> > readDATAtoVectors(string totalFname, int N); //moved to beginning of the file

vector<double> Ising::mean_S_corrrelation(string totalFname) {
    /* Domain : calculation of physical quantities.
     * Calculate macroscopic connected spin-spin correlation function obtained
     *  by taking time mean value of instantaneous correlations.
     */

    //Raw simulation data flows from disk to vector container
    vector<vector<short> > SS = readDATAtoVectors(totalFname, N);
    
    //Declaration of iterators
    vector<vector<short> >::iterator t;
    vector<short> s;
    vector<short>::iterator pos;
    
    //Local container for instant correlation
    vector<double> S_corr(maxGdist, 0);
    
    cout<<"\n\n    MEASURE OF THE SPIN-SPIN CORRELATION FUNCTION"<<endl<<endl;
    cout << "\n Number of measures = " << SS.size() << endl; //Size of SS should be equal to number of measures (check)
    
    int i;

    for (t = SS.begin(); t != SS.end(); t++) { //loop over time "t" (time unity in MC cycles = prod_t/measure_f)

        i = 0;
        for (pos = t->begin(); pos != t->end(); pos++) { //loop over chain sites "pos" at fixed time
            S[i] = *pos;
            i++;
        } //We are leaving the loop with an instantaneous state of the chain contained in the Ising class member vector S

        double mcycle = order_parameter();
        S_corr = S_correlation(maxGdist, mcycle); 
        for (int j = 0; j < maxGdist; j++)
            Gs[j] += S_corr[j]; //accumulation of instantaneous results (we are in time loop)
    }

    for (int j = 0; j < maxGdist; j++)
        Gs[j] /= SS.size(); //normalization by number of measures

    //DISPLAY SECTION (from here to the end of the function)
    cout << "\nFinal order parameter <S> = " << order_parameter();
    cout << "\n\nSpin-spin correlation function";
    for (int j = 0; j < maxGdist; j++)
        if (Gs[j] > 0)
            cout << "\nG[" << j << "]=" << Gs[j];

    cout << "\n\nTemperature calculated from simulations (*)";

    for (int j = 1; j < maxGdist; j++)
        if (Gs[j] > 0)
            cout << "\nT[" << j << "]=" << 1 / atanh(Gs[j] / Gs[j - 1]);
    
    cout << "\n--------------------";
    cout << "\n(*) We use the fact that correlation is equal to [th(1/T)]^dist";
    cout << "\nwhere \"dist\" represents distance between spins)";
    cout << "\nValues above corresponds to 1/atanh(G[dist] / G[dist-1]), ";
    cout << "\nThis gives the measured value of the temperature T";
    
    return Gs;
}

