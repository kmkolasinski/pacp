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

//#include<stack> //LIFO (Last In First Out) container for Wolff MC algorithm

#include "Ising5_2.h"
//#include "Tools.h"


using namespace std;

string gettotalFname(); //forward (Tools)
void displayDATAfromVectors(string totalFname, int N); //forward (Tools)

int main(int argc, char** argv) {

    //SIMULATION INPUT DATA
    int chainLength = 2000;
    double temperature = 0.8;
    double magneticField = 0.0;
    int therm_t = 500; //thermalization time
    int measure_f = 5; //measure frequency
    int prod_t = 10000; //productiontime
    int maxGdist = 10; //domain of correlation function extends from 0 to maxGdist-1
    int maxTsep = 2;
    double initState = 0.7; //values between 0.5 and 1.0
    //0.5 --> high temperature limit --> totally disordered state, <S>=0.
    //1.0 --> low temperature limit --> totally ordered state, <S>=1.
    //***When needed add other initializations here

    //CREATION AND INITIALIZATION OF THE INSTANCE OF ISING CLASS
    Ising chain(chainLength, temperature, magneticField, initState, maxGdist, maxTsep);


    //CALL TO THE ISING CLASS MEMBER FUNCTION PERFORMING SIMULATION
    chain.MC_simulation(therm_t, prod_t, measure_f);

    
    //Exemple from section D of implementation file (Ising5_2.cpp)
    chain.mean_S_corrrelation(gettotalFname().c_str());
    
    //Check of data. To be used for small time x space values
    //displayDATAfromVectors(gettotalFname().c_str(),chainLength); 
    
    cout << "\n\n     T H E   E N D " << endl;

    return 0;
}

