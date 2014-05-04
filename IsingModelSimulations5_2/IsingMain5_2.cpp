/* 
 * File:   IsingMain5_2.cpp
 *
 */

//////////////////////////////////////////////////////////////////////////////
// MAIN FUNCTION
//////////////////////////////////////////////////////////////////////////////

#include <stdio.h>      
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>		
#include <cstring>		
#include <iomanip>
#include <cmath>

//#include<stack> //LIFO (Last In First Out) container for Wolff MC algorithm

#include "Tools.h"
#include "Ising5_2.h"
#include "Ising2D.h"
#include "Tests/IsingTestChi.h"
#include "Tests/IsingTestCC2D.h"



using namespace std;



int main(int argc, char** argv) {

    
    // IsingTestChi isingTestChi;   //  Ising1D test of Chi in function of T
    // IsingTestCC2D isingTestCC2D; //  Ising2D test of CC in function of T for different lattices
    
    
    //SIMULATION INPUT DATA
    int chainLength = 500;
    double temperature = 1.0;
    double magneticField = 0.0;
    int therm_t = 500; //thermalization time
    int measure_f = 10; //measure frequency
    int prod_t = 100000; //productiontime
    int maxGdist = 10; //domain of correlation function extends from 0 to maxGdist-1
    int maxTsep = 2;
    double initState = 0.7; //values between 0.5 and 1.0
    //0.5 --> high temperature limit --> totally disordered state, <S>=0.
    //1.0 --> low temperature limit --> totally ordered state, <S>=1.
    //***When needed add other initializations here    
    
    
    if(false){
    // ------------------------------------------------------------------------
    // test w funkcji production time
    // ------------------------------------------------------------------------
    ofstream data_out("error_test.txt");
    Ising chain(chainLength, temperature, magneticField, initState, maxGdist, maxTsep);
    for ( int prod_t = 1000 ; prod_t < 100000 ; prod_t *= 2  ){
        chain.MC_simulation(therm_t, prod_t, measure_f);
        double error=chain.ERROR(gettotalFname().c_str(), ERROR_CHI);
        data_out << std::scientific <<  prod_t << "\t" << error  << "\n";
        cout   << std::scientific <<  prod_t << "\t" << error  << "\n";
    }
    data_out.close();
    }
    
    cout << "\n\n     T H E   E N D " << endl;

    return 0;
}

