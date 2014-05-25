

#ifndef ISINGTESTLATTICEEVOLUTION_H
#define	ISINGTESTLATTICEEVOLUTION_H


#include "IsingTest.h"
#include <string>

class IsingTestLatticeEvolution : public IsingTest {

    public:
        IsingTestLatticeEvolution():IsingTest(){
            test_name = "Observation of lattice evolution";   
            test_info = string(" Run test to look for lattice evolution through temperature.  \n");            
            info(); 
            run(); // run all calculations
            //string cmd =string("cd ")+test_dir_output+string(";gnuplot IsingTestLatticeEvolution.plt");
            //int info = system(cmd.c_str());

        }   
        
        /**
         * Search spin lattice for clusters.
         */
        void run(){
         
        cout << " Running test..." << endl;                
        cout << " God, shave the Queen..." << endl;
        
        //SIMULATION INPUT DATA
        int chainLength      = 100;    //initial lattice size
        double magneticField = 0.0;
        double T             = 5.0;
        int therm_t          =  50;    //thermalization time
        int measure_f        =   1;    //measure frequency
        int prod_t           =1000;    //productiontime
        int maxGdist         =  10;    //domain of correlation function extends from 0 to maxGdist-1
        int maxTsep          =   2;
        double initState     = 0.7;    //values between 0.5 and 1.0
 
        

        //int info = system(cmd.c_str());
        int info;       
        Ising2D lattice(chainLength, T, magneticField, initState, maxGdist, maxTsep);        
        lattice.MC_simulation(therm_t, prod_t, measure_f,METHOD_WOLFF);
        char temp;
        for(int k=0 ; k<30 ; k++)  {
            //string s = static_cast<ostringstream*>( &(ostringstream() << k) )->str();            
            temp=k+65;
            string fileout = test_dir_output+"lat/IsingTestLatticeEvolution"+temp+".txt";
            ofstream data_out(fileout.c_str()); 
            for(int i=0 ; i<chainLength ; i++){
                for(int j=0; j< chainLength; j++) data_out << lattice.S[j+chainLength*i] << "\t";
                data_out << "\n";
            }
            cout << "zapis\n" << endl;
            data_out.close();
            lattice.T-=0.1;
            for(int i=0 ; i<therm_t ; i++) lattice.cycle();
        }
        
                        
        }// end of run()
        

};


#endif	/* ISINGTESTLATTICEEVOLUTION_H */

