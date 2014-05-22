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




#include "Ising5_2.h"
#include "Tools.h"



using namespace std;

long int Ising::seed1 = 46723;
long int Ising::seed2 = 2593459;

vector<vector<short> > readDATAtoVectors(string totalFname, int N); //forward declaration (Tools)




//////////////////////////////////////////////////////////////////////////////
//  section A             CONSTRUCTOR OF OBJECT "ISING"                     //          
//////////////////////////////////////////////////////////////////////////////

Ising::Ising(int NN, double TT, double hh, double p, int maxGd, int maxTsep) {
    //CONSTRUCTOR FOR 1 DIMENSIONAL ISING CHAIN
    N = NN; //set length of the chain
    T = TT; //set temperature
    h = hh; //set magnetic field
    maxGdist          = maxGd;
    maxTimeseparation = maxTsep;
    ising_method_type      = METHOD_METROPOLIS; // default method for calulation

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

vector<short> Ising::Wolff_cycle() {
    /* Domain : Monte Carlo simulation
     * Perform N Wolff updates of cluster of spins 
     * i.e. one MC cycle or sweep*/
    stack<int> border;
    int oldspin, newspin;
    double padd = 1 - exp(-2 / T); //probability of adding a site to border
    int site, siteL, siteR;


    for (int i = 0; i < N; i++) {

        site = (int) (N * RNG::drandom(seed2));

        border.push(site);
        oldspin = S[site];
        newspin = -S[site];
        S[site] = newspin;

        while (!border.empty()) {
            site = border.top();
            border.pop();

            siteL = (site - 1 + N) % N;
            if (S[siteL] == oldspin)
                if (RNG::drandom(seed2) < padd) {
                    border.push(siteL);
                    S[siteL] = newspin;
                }

            siteR = (site + 1) % N;
            if (S[siteR] == oldspin)
                if (RNG::drandom(seed2) < padd) {
                    border.push(siteR);
                    S[siteR] = newspin;
                }

        }
    }
    return S;
}


vector<short> Ising::cycle(){

    switch(ising_method_type){
        case  METHOD_METROPOLIS: return Metropolis_cycle();break;
        case  METHOD_WOLFF:      return Wolff_cycle();break;
        default: return Metropolis_cycle();
    }
    
}


void Ising::MC_simulation(int therm_t, int prod_t, int measure_f,ISING_METHOD_TYPE mt) {
    /* Domain : Metropolis simulation of Ising 1d model
     * Perform desired number of MC cycles
     * therm_t = thermalization cycles
     * prod_t  = production cycles
     * Store raw results (states of the chain) in the disk file. */

    ising_method_type = mt;    
    RNG::initializeRNG(seed1); // initializes integer random number generator

    //*** INITIAL MC CYCLES ***
    //Thermalization of initial non-equilibrium state.
    for (int t = 0; t < therm_t; t++) {
        cycle();
    }

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
    for (int t = 0; t < prod_t; t++) {
        if (!(t % measure_f)) {
            s = cycle();
            for (int i = 0; i < N; i++)
                DATA << s[i] << " ";
            DATA << endl;
        } else
            cycle();
        /*
         Please make copy of this file (Ising5_2.cpp)
         */
        
    }
    DATA.close();
}





//////////////////////////////////////////////////////////////////////////////
// section C                  DATA ANALYSIS                                 //          
//////////////////////////////////////////////////////////////////////////////

double Ising::E() {
    double E;
    E = 0;
    for (int i = 0; i < N; i++) {
        if (i == (N - 1)) {
            E += -S[i] * S[0];
        } else {
            E += -S[i] * S[i + 1];
        }
    }
    return E / (N);
}

vector<double> Ising::E_correlation(int maxGdist, double E) {
    /*Domain: data analysis, calculation of physical quantities
     * Calculate instantaneous connected correlation function
     * E represents instantenous (i.e. lattice) mean value of energy. 
     *         Ge(x)= <E(0)E(x)> - <E>^2
     */

    double E2 = E*E;
    vector<double> corr(maxGdist,0);
    for (int dist = 0; dist < maxGdist; dist++) {
        for (int pos = 0; pos < N; pos++) {
            corr[dist] += S[pos] * S[(pos + dist) % N]-E2;
        }
        corr[dist] /= N;
    }
    return corr;

}

vector<double> Ising::mean_E_correlation(vector <vector<short> > SS) {
    /*
     * Function allows to measure energy correlation by calculating its
     * time mean value.
     * 
     */
     
    //Declaration of iterators
    vector<vector<short> >::iterator t;
    vector<short> s;
    vector<short>::iterator pos;

    //TEMPORARY VARIABLE FOR E CORRELATION
    vector<double> E_corr(maxGdist,0);

    int i;

    for (t = SS.begin(); t != SS.end(); t++) { //loop over time "t" (time unity in MC cycles = prod_t/measure_f)

        i = 0;
        for (pos = t->begin(); pos != t->end(); pos++) { //loop over chain sites "pos" at fixed time
            S[i] = *pos;
            i++;
        } //We are leaving the loop with an instantaneous state of the chain contained in the Ising class member vector S

        double E_tot =E();
        E_corr = E_correlation(maxGdist, E_tot);
        for (int j = 0; j < maxGdist; j++)
            Ge[j] += E_corr[j]; //accumulation of instantaneous results (we are in time loop)
    }
    for (int j = 0; j < maxGdist; j++)
        Ge[j] /= SS.size(); //normalization by number of measures

    return Ge;
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
    double meanE, meanEsq;

    meanE = 0;
    meanEsq = 0;
    int n=0;
    for (int t = 0; t < 100000; t++) {        
     
            cycle();
            double ee = E();            
            meanE  += ee;
            meanEsq+= ee*ee;
            n++;
    }
    meanE /= n;
    meanE *= meanE;
    meanEsq /= n;

    return (meanEsq - meanE) / T / T / N; //per site
}

double Ising::CC(string totalFname) {
    // Calculates the specific heat 
    // which is stored in file totalFname
    // parameter: totalFname - name of file    
    //Raw simulation data flows from disk to vector container
    vector<vector<short> > SS = readDATAtoVectors(totalFname, N);
    return CC(SS);
}

double Ising::CC(vector<vector<short> > SS) {
    /*
     * returns value of the specific heat according to relation:
     * C=bheta/T <E^2>-<E>^2
     * SS - is the raw data vector stored in memory
     * this function needs E() to work properly.
     */
    double meanE, meanEsq;
    int t_steps = SS.size();

    meanE = 0;
    meanEsq = 0;
    for (int j = 0; j < t_steps; j++) {
        for (int i = 0; i < N; i++) {
            S[i] = SS[j][i];
        }
        double ee = E();
        meanE += ee;
        meanEsq += ee*ee;
    }

    meanE /= t_steps;
    meanE *= meanE;
    meanEsq /= t_steps;

    return (meanEsq - meanE) / T / T / N;
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
    aM2=aM2/t_steps/N/N;
    
    return (aM2-aM*aM)/T; //Chi
}


double Ising::aveE(string totalFname) {
    // Calculates the average value of Energy over all simulation time

    //Raw simulation data flows from disk to vector container
    vector<vector<short> > SS = readDATAtoVectors(totalFname, N);
    //number of time samples
    int t_steps = SS.size();
    //auxiliary variables
    double mean_E = 0;
    for (int j = 0; j < t_steps; j++) {
        for (int i = 0; i < N; i++) {
            S[i] = SS[j][i];
        }
        double ee = E();
        mean_E += ee;
    }
    return mean_E / t_steps;
}

cluster_stat Ising::cluster_freq1D(string totalFname) {
    // Calculation of the cluster size distribution.
    // Input: name of the file containing the time evolution of the system created by MC simulation.
    // Output: cluster_stat structure     //parameter: totalFname - name of file.
    // This function feed in data mean_cluster_freq1D and forward its result to output.

    //Raw simulation data flows from disk to vector container SS
    vector<vector<short> > SS = readDATAtoVectors(totalFname, N);

    return mean_cluster_freq1D(SS);
}
cluster_stat Ising::mean_cluster_freq1D(string totalFname) {
    // Calculation of the cluster size distribution
    // Input: time evolution of the system created by MC simulation
    // i.e. SS represents raw simulation data
    // Output: cluster_stat structure 
    vector<vector<short> > SS = readDATAtoVectors(totalFname, N);
    return mean_cluster_freq1D(SS);
}

cluster_stat Ising::mean_cluster_freq1D(vector<vector<short> > SS) {
    // Calculation of the cluster size distribution
    // Input: time evolution of the system created by MC simulation
    // i.e. SS represents raw simulation data
    // Output: cluster_stat structure 
 
    cluster_stat instantClstat;
    for (int j = 0; j < N + 1; j++) CFD.CF.push_back(0);
    CFD.cmax = 0;

    //number of time samples
    int t_steps = SS.size();
    for (int j = 0; j < t_steps; j++) {
        instantClstat = cluster_freq1D(SS[j]);
        for (int i = 1; i < N + 1; i++)
            CFD.CF[i] += instantClstat.CF[i];
        if (instantClstat.cmax > CFD.cmax)
            CFD.cmax = instantClstat.cmax;
    }
    //normalization by number of measures and sites
    //additional multiplying by factor i - size of each normalized cluster number
    //gives us vector structure of probabilities of observing clusters of given size
double temp = t_steps * N;
for (int i = 1; i < N + 1 ; i++)
    CFD.CF[i] *= i/temp;

return CFD;
}

cluster_stat  Ising::cluster_freq1D(vector<short> SS) {
    /*Domain : data analysis: calculation of the cluster size distribution
     * For 1-dimensional system the T=0 corresponds to criticality
     * When approaching critical state the system develops large clusters of 
     * spins pointing in the same directions. Such clusters are named domains.
     * When going closer and closer to critical temperature (T=0 in 1d case)
     * clusters are not only growing in size but their distribution flatten
     * meaning that domains of any size are present close to criticality.
     * The function below is designed to check this behavior. It measures
     * the distribution of domain for instant chain state. Function returns structure
     * cluster_stat containing two members: 
     * - vector CF: CF[L] contains the frequency of domain of size L
     * - cmax : contains the size of the largest domain in the distribution
     * This function returns the instantaneous cluster distribution
     */
    vector<float> tt(N + 1, 0); //tt is created to initialize the vector component
    //of the cls instance of the cluster_stat structure below
    cluster_stat cls = {tt, 1}; //container for results
    int pos = 0; //initial current position on the chain
    int sg = SS[pos]; //sign of the spin at current position (here of the spin at the beginning of the chain)
    int cl_L = 0; //length of the current cluster (initialized to zero)

    //If the first and the last spin of the chain point in the same direction
    //they belong to the same cluster for the chain with periodic boundary conditions
    //the block below and the last block treat such situations

    int first_cl_L;
    if (SS[0] == SS[N - 1]) { //Special case: cluster wrap around the border.
        while (SS[pos] == sg) {
            pos++;
            cl_L++;
        }
        first_cl_L = cl_L;
        sg = -sg;
        if (cls.cmax < cl_L) cls.cmax = cl_L; //I am bigger then you, ... no!, I am bigger then you.
    }

    if (cls.cmax == N) //The special case occurs when the whole chain is in unique domain state.
        cls.CF[N] = 1; //In this case the treatment is accomplished and we return result.
    else { //Here the most common treatment is implemented. 
        while (pos < N) { // We progress till the end of the chain,
            cl_L = 0;
            while (SS[pos] == sg) { //counting the size of one domain after another.
                pos++;
                cl_L++;
            }
            cls.CF[cl_L]++;
            sg = -sg;
            if (cls.cmax < cl_L) cls.cmax = cl_L;
        }

        int last_cl_L = cl_L; //However we keep in mind that the last domain 
        if (SS[0] == SS[N - 1]) { //may glue to the first one due to the periodic boundary condition.
            cls.CF[last_cl_L]--; //In this case the last domain should not be counted
            cls.CF[first_cl_L + last_cl_L]++; // because the length of crossing boundary domain is the joint length
            if (cls.cmax < first_cl_L + last_cl_L) //of the first and the last "domains".
                cls.cmax = first_cl_L + last_cl_L;
        }
    }
    return cls;
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

    int i = 0;

    for (t = SS.begin(); t != SS.end(); t++) { //loop over time "t" (time unity in MC cycles = prod_t/measure_f)
        for (pos = t->begin(); pos != t->end(); pos++) { //loop over chain sites "pos" at fixed time
            mt[i] += (double) (*pos);
        }
        mt[i] /= (double) N;
        i++;
        //if(i==maxTimeseparation) break;
    }


    cout << "\n--------------------";

    for (int j = 0; j < maxTimeseparation; j++)
        cout << "\nmt[" << j << "]=" << mt[j];

    cout << "\n--------------------";



    double meanm = 0.;
    for (int j = 0; j < maxTime; j++)
        meanm += mt[j];
    meanm /= (double) maxTime;


    /////////////////// MAIN LOOP ///////////////////////////////////
    for (int j = 0; j < maxTimeseparation; j++) {

        for (int k = 0; k < maxTime; k++) {

            int index = j + k;
            if ((j + k) >= maxTime) index -= maxTime;

            Gt[j] += (mt[k] * mt[index] - meanm * meanm);

        }

    }
    /////////////////////////////////////////////////////////////////

    double Gt_int = 0.;
    fstream DATA("Gt.dat", ios::out);
    for (int j = 0; j < maxTimeseparation; j++) {
        cout << "\nGt[" << j << "]=" << Gt[j];
        DATA << Gt[j] << "\n";
        Gt_int += Gt[j];
    }
    cout << "\n--------------------";

    cout << "i = " << i << endl;
    cout << "meanm = " << meanm << endl;
    cout << "Gt integral = " << Gt_int << endl;
    Gt_int /= Gt[0];
    cout << "tau = " << Gt_int << endl;
    Gt_int = 1. / Gt_int;
    cout << "1/tau = " << Gt_int << endl;

    return Gt;

}

//function returning sigma for specific error_type (see enum file in *.h file)

double Ising::ERROR(string totalFname, ISING_ERROR_TYPE error_type) {
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

    double observable;
    double sum2 = 0, sum1 = 0;
    int m = 100;

    for (int j = 0; j < m; j++) {

        // we choose random set (with returns)
        for (int i = 0; i < n; i++) {
            int time_step = (int) (n * RNG::drandom(seed2));
            SS[i] = SSdata[time_step];
        }
        switch (error_type) { // choose of which observable error is calculated
            case ERROR_CHI:
                observable = Chi(SS);
                break;
            case ERROR_CC:
                observable = CC(SS);
                break;
            case ERROR_OP:
                observable = order_parameter(SS);
                break;
        }

        sum2 += observable*observable;
        sum1 += observable;
    }
    sum2 /= m;
    sum1 /= m;

    double sigma = sqrt(sum2 - sum1 * sum1);

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
    int pos = 0; //initial current position in the chain
    int sg = S[pos]; //initial sign of current position (i.e. spin orientation)
    int cl_L = 1; //initial length of the current cluster
    int max = 0; //temporary variable used in the for loop
    pos = 0;
    bool nbound_c = 1;
    if (sg == S[N - 1]) nbound_c = 0;
    int temp = 0; //first domain check
    int sg_size = 0; //size of first domain
    while (pos < N) {
        if (nbound_c) {
            if (pos != (N - 1)) {
                if (S[pos] == S[++pos])cl_L++;
                else {
                    tt[cl_L - 1] += 1;
                    if (max < cl_L) max = cl_L;
                    cl_L = 1;
                }
            } else {
                tt[cl_L - 1] += 1;
                if (max < cl_L) max = cl_L;
                cl_L = 1;
                pos++;
            }
        } else {
            if (pos != (N - 1)) {
                if (S[pos] == S[++pos])cl_L++;
                else {
                    temp++;
                    if (temp == 1) sg_size = cl_L;
                    tt[cl_L - 1] += 1;
                    if (max < cl_L) max = cl_L;
                    cl_L = 1;
                }
            } else {
                if (cl_L != N) {
                    cl_L += sg_size;
                    tt[sg_size - 1] -= 1;
                }
                tt[cl_L - 1] += 1;
                if (max < cl_L) max = cl_L;
                cl_L = 1;
                pos++;
            }
        }
    }
    cls.CF = tt;
    cls.cmax = max;

    return cls;
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

double Ising::order_parameter(string totalFname) {
    // Calculates the magnetization (order parameter) of the system
    // which is stored in file totalFname
    //parameter: totalFname - name of file

    //Raw simulation data flows from disk to vector container
    vector<vector<short> > SS = readDATAtoVectors(totalFname, N);

    return order_parameter(SS);

}

double Ising::order_parameter(vector<vector<short> > SS) {
    // Calculates the magnetization (order parameter) of the system
    // which is stored in vector SS
    //parameter: SS - Raw simulation data

    //number of time samples
    int t_steps = SS.size();
    //auxiliary variables
    double M_temp, aM;

    aM = 0.;
    for (int j = 0; j < t_steps; j++) {
        M_temp = 0.;
        //calculation of magnetisation M
        for (int i = 0; i < N; i++) {
            M_temp += SS[j][i];
        }
        //calculation of <M>
        aM += M_temp;
    }
    //normalization by number of measures and sites
    aM = aM / t_steps / N;

    return aM; //M


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

vector<double> Ising::mean_S_correlation(string totalFname) {
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

    cout << "\n\n    MEASURE OF THE SPIN-SPIN CORRELATION FUNCTION" << endl << endl;
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

vector<double> Ising::mean_S_correlation(vector <vector<short> > SS) {
    /* Domain : calculation of physical quantities.
     * Calculate macroscopic connected spin-spin correlation function obtained
     *  by taking time mean value of instantaneous correlations.
     */

    //Declaration of iterators
    vector<vector<short> >::iterator t;
    vector<short> s;
    vector<short>::iterator pos;
;

    //Local container for instant correlation
    vector<double> S_corr(maxGdist, 0);

    cout << "\n\n    MEASURE OF THE SPIN-SPIN CORRELATION FUNCTION" << endl << endl;
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

    /*
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
    */
    return Gs;
}

vector<short> Ising::Domain_snapshot(){
    vector<short> domains(S);
    int oldspin;
    int site=0;
    int esd,osd;//even and odd size of domain
    esd=2;//corresponds to +1 spins
    osd=3;//corresponds to -1 spins        
    for (int i = 0; i < N; i++) {
        oldspin = domains[site];
        for(int j = 0; j < N; j++){
            if(domains[i]==oldspin){
                oldspin==1? domains[i]=esd : domains[i]=osd;
            }else{
                
            }
        }
    }
    return S; 
};

