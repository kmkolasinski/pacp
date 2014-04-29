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

#include "Ising5_2.h"
//#include "Tools.h"


bool bTestChiTemp       = false; // test caclulate Chi(T)
bool bErrorProdTimeTest = false; // test caclulate Errro of Chi in function of prod_t

using namespace std;

string gettotalFname(); //forward (Tools)
void displayDATAfromVectors(string totalFname, int N); //forward (Tools)

int main(int argc, char** argv) {

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
    
    if(bTestChiTemp){
        // ------------------------------------------------------------------------
        // Chi in function of Temperature
        // ------------------------------------------------------------------------


        ofstream data_out("Temp_test.txt");
        for( double temperature = 0.3 ; temperature < 3.0 ; temperature += 0.05  ){
        //CREATION AND INITIALIZATION OF THE INSTANCE OF ISING CLASS
        Ising chain(chainLength, temperature, magneticField, initState, maxGdist, maxTsep);


        //CALL TO THE ISING CLASS MEMBER FUNCTION PERFORMING SIMULATION
        chain.MC_simulation(therm_t, prod_t, measure_f);


        //Exemple from section D of implementation file (Ising5_2.cpp)
        //chain.mean_S_corrrelation( gettotalFname().c_str());
        //cout<< endl << chain.Chi(gettotalFname().c_str()) << " " << exp(2./temperature)/temperature << endl;

        double chi = chain.Chi(gettotalFname().c_str());
        double chi_anal = exp(2./temperature)/temperature;
        double error=chain.ERROR( gettotalFname().c_str() );

        data_out << std::scientific <<  temperature  << "\t"
                                         << chi_anal << "\t"
                                         << chi      << "\t"
                                         << error    << "\t"
                                        << chi-error << "\t"
                                        << chi+error << endl;

        cout << std::scientific <<  temperature  << "\t" << chi << "\t" << chi_anal << "\t" <<  error << endl;

        }
        data_out.close();
    } // end of if(bTestChiTemp)
    
    if(bErrorProdTimeTest){
    // ------------------------------------------------------------------------
    // test w funkcji production time
    // ------------------------------------------------------------------------
    ofstream data_out("error_test.txt");
    Ising chain(chainLength, temperature, magneticField, initState, maxGdist, maxTsep);
    for ( int prod_t = 1000 ; prod_t < 100000 ; prod_t *= 2  ){
        chain.MC_simulation(therm_t, prod_t, measure_f);
        double error=chain.ERROR( gettotalFname().c_str() );
        data_out << std::scientific <<  prod_t << "\t" << error  << "\n";
        cout   << std::scientific <<  prod_t << "\t" << error  << "\n";
    }
    data_out.close();
    }
    
    
    //Check of data. To be used for small time x space values
    //displayDATAfromVectors(gettotalFname().c_str(),chainLength); 
    
    cout << "\n\n     T H E   E N D " << endl;

    return 0;
}

