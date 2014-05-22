

#ifndef ISINGTESTCLUSTERFREQ2D_H
#define	ISINGTESTCLUSTERFREQ2D_H


#include "IsingTest.h"

class IsingTestClusterFreq2D : public IsingTest {

    public:
        IsingTestClusterFreq2D():IsingTest(){
            test_name = "Observation of cluster size distribution";   
            test_info = string(" Run test to look for cluster size distribution.  \n");
            
            info(); 
            run(); // run all calculations
            string cmd =string("cd ")+test_dir_output+string(";gnuplot IsingTestClusterFreq2D.plt");
           // int info = system(cmd.c_str());
        }   
        
        /**
         * Search spin lattice for clusters.
         */
        void run(){
         
        cout << " Running test..." << endl;                
        cout << " God, shave the Queen..." << endl;
        
        //SIMULATION INPUT DATA
        int chainLength      = 50;    //initial lattice size
        double magneticField = 0.0;
        double T             = 10.0;
        int therm_t          =  50;    //thermalization time
        int measure_f        =   1;    //measure frequency
        int prod_t           =1000;    //productiontime
        int maxGdist         =  10;    //domain of correlation function extends from 0 to maxGdist-1
        int maxTsep          =   2;
        double initState     = 0.7;    //values between 0.5 and 1.0
 
        string fileout = test_dir_output+"IsingTestClusterFreq2D.txt";
        ofstream data_out(fileout.c_str()); 
                
        Ising2D lattice(chainLength, T, magneticField, initState, maxGdist, maxTsep);        
        lattice.MC_simulation(therm_t, prod_t, measure_f,METHOD_WOLFF);
        
        /*for(int i=0 ; i<chainLength ; i++){
            for(int j=0; j< chainLength; j++){
                data_out << lattice.S[j+chainLength*i] << "\t";
            }
            data_out << "\n";
        }*/
        vector <short> cf(lattice.cluster_freq2D().SS);        
        for(int i=0 ; i<chainLength ; i++){
            for(int j=0; j< chainLength; j++) data_out << cf[j+chainLength*i] << "\t";
            data_out << "\n";
        }
        
        
        data_out.close();
                
        }// end of run()

};


#endif	/* ISINGTESTCLUSTERFREQ2D_H */

