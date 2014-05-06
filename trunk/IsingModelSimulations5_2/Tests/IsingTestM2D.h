/* 
 * File:   IsingTestM2D.h
 * Author: edyta
 *
 * Created on 5 maj 2014, 15:22
 */

#ifndef ISINGTESTM2D_H
#define	ISINGTESTM2D_H


#include "IsingTest.h"

class IsingTestM2D : public IsingTest {

    public:
        IsingTestM2D():IsingTest(){
            test_name = "M in function of T .";   
            test_info = string(" Run test to calculate the M value in function of  \n")+
                        string(" temperature T for the 2D Ising lattice. The test is performed \n")+
                        string(" for different lattice sizes: 2x2, 4x4, 8x8, 16x16. Critical \n")+
                        string(" temperature is Tc=2.27.\n")+
                        string(" The results are saved in proper files:")+  
                        string(" 'IsingTestM2D.txt' and 'IsingTestM2D.png' in 'tests_out' directory. \n")+
                        string(" See run() function for more details.");
            
            info(); 
            run(); // run all calculations
            string cmd =string("cd ")+test_dir_output+string(";gnuplot IsingTestM2D.plt");
            int info = system(cmd.c_str());
        }   
        
        /**
         * Tests Chi in function of Temperature for 1D model
         */
        void run(){
         
        cout << " Running test..." << endl;                
        cout << " It may take some time..." << endl;
        //SIMULATION INPUT DATA
        int chainLength      = 2;        
        double magneticField = 0.0;
        int therm_t          = 1000;   //thermalization time
        int measure_f        = 1;     //measure frequency
        int prod_t           =10000;  //productiontime
        int maxGdist         = 10;     //domain of correlation function extends from 0 to maxGdist-1
        int maxTsep          = 2;
        double initState     = 0.7;    //values between 0.5 and 1.0
    
                
        vector<double> dataNxN[4];
        vector<double> dataT;
        for(double T = 1.0 ; T < 5.0 ; T += 0.1 ){ 
            dataT.push_back(T);
        }
        // ---------------------------------------------------------
        // Calculation for NxN lattice.
        // ---------------------------------------------------------
        for(int i=0 ; i < 4 ; i++){                    
            chainLength  = pow(2,i+1); // choosing proper lattice size
            cout << "Starting simulation for lattice "<<chainLength << "x" <<  chainLength << "...\n";  
            
            for(int t = 0 ; t < dataT.size() ; t++ ){ 
                double T = dataT[t];
                cout << "T=" << T << "\t:";
                Ising2D lattice2d(chainLength, T, magneticField, initState, maxGdist, maxTsep);        
                lattice2d.MC_simulation(therm_t, prod_t, measure_f,METHOD_METROPOLIS);
                double M = lattice2d.order_parameter(gettotalFname().c_str());
                dataNxN[i].push_back(M);
            }// end of for            
        }
        
        // saving data to file
        string fileout = test_dir_output+"IsingTestM2D.txt";
        ofstream data_out(fileout.c_str());     
        
        for(int t = 0 ; t < dataT.size() ; t++ ){
             data_out << std::scientific << dataT[t] << "\t";   
             for(int i = 0; i < 4 ; i++){
                 data_out << fabs(dataNxN[i][t]) << "\t";
             }
             data_out << endl;
        }
        
        
        }// end of run()

};



#endif	/* ISINGTESTM2D_H */

