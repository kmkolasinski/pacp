/*
* File:   IsingTestTCorr.h
* Author: mdyndal
*
*/

#ifndef ISINGTESTTCORR_H
#define ISINGTESTTCORR_H

#include "IsingTest.h"

class IsingTestTCorr : public IsingTest {

   public:
       IsingTestTCorr():IsingTest(){
           test_name = "Test of the autocorrelation function Psi";  
           test_info = string(" Run test to calculate the autocorrelation function Psi as a function of time");
           
           info();
           run(); // run all calculations
       }  
       

       void run(){
       
       cout << " Running test..." << endl;    

        //SIMULATION INPUT DATA
        int chainLength      = 500;        
        double magneticField = 0.0;
        int therm_t          = 1000;   //thermalization time
        int measure_f        = 10;     //measure frequency
        int prod_t           = 10000;  //productiontime
        int maxGdist         = 10;     //domain of correlation function extends from 0 to maxGdist-1
        int maxTsep          = 2;
        double initState     = 0.7;    //values between 0.5 and 1.0


        cout << " Start  T=" << 0.3 << endl;
        cout << " End    T=" << 2.0 << endl;
        cout << " Step  dT=" << 0.1 << endl;
        cout << " It may take some time..." << endl;
        cout << "----------------------------------------------------" << endl;
        
        for( double temperature = 0.3 ; temperature < 2.0 ; temperature += 0.1  ){        
        //CREATION AND INITIALIZATION OF THE INSTANCE OF ISING CLASS
        Ising chain(chainLength, temperature, magneticField, initState, maxGdist, maxTsep);
        //CALL TO THE ISING CLASS MEMBER FUNCTION PERFORMING SIMULATION
        chain.MC_simulation(therm_t, prod_t, measure_f);
        chain.t_correlation(prod_t/measure_f, gettotalFname().c_str());
        }// end of for(temp)
                                       
        }// end of run()


};


#endif  /* ISINGTESTTCORRH */

