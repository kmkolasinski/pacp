
#ifndef ISINGTESTSS_CORR_H
#define	ISINGTESTSS_CORR_H

#include "IsingTest.h"

class IsingTestSS_cor : public IsingTest {
    
    public:
        IsingTestSS_corr():IsingTest(){
            test_name = "Spin-spin correlation test.";   
            test_info = string(" Run test to calculate the Chi value in function of \n")+
                        string(" temperature T for the 1D Ising chain.\n The results are saved in proper files:")+  
                        string(" 'IsingTestChi1D.txt' and 'IsingTestChi1D.png' in 'tests_out' directory. \n")+
                        string(" See run() function for more details.");
            
            info(); 
            run(); // run all calculations
            string cmd =string("cd ")+test_dir_output+string(";gnuplot IsingTestChi1D.plt");
            int info = system(cmd.c_str());
        }   
        void run(){
         
        cout << " Running test..." << endl;    
            
        //SIMULATION INPUT DATA
        int chainLength      = 5000;        
        double magneticField = 0.0;
        int therm_t          = 1000;   //thermalization time
        int measure_f        = 10;     //measure frequency
        int prod_t           = 10000;  //productiontime
        int maxGdist         = 10;     //domain of correlation function extends from 0 to maxGdist-1
        int maxTsep          = 2;
        double initState     = 0.7;//values between 0.5 and 1.0
        double minT          = 0.4;
        double maxT          = 1.5;
        
        string fileout = test_dir_output+"IsingTestSS_corr.txt";
        ofstream data_out(fileout.c_str());
        
        for( double temperature = minT ; temperature < maxT ; temperature += 0.1  ){         
            //CREATION AND INITIALIZATION OF THE INSTANCE OF ISING CLASS
            Ising chain(chainLength, temperature, magneticField, initState, maxGdist, maxTsep);
            //CALL TO THE ISING CLASS MEMBER FUNCTION PERFORMING SIMULATION
            chain.MC_simulation(therm_t, prod_t, measure_f);
            double op=chain.order_parameter();
            vector <double> SS_corr = chain.S_correlation(maxGdist,op);
            for(int i=0;i<SS_corr.size;i++){
                data_out<< i << ";" << SS_corr[i] << endl;
            }
        }
        data_out.close();                                    
        }

};


#endif	/* ISINGTESTCHI_H */

