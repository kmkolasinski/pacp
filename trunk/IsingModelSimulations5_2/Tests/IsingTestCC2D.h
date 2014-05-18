

#ifndef ISINGTESTCC2D_H
#define	ISINGTESTCC2D_H


#include "IsingTest.h"

class IsingTestCC2D : public IsingTest {

    public:
        IsingTestCC2D():IsingTest(){
            test_name = "Test of specific heat in function of T .";   
            test_info = string(" Run test to calculate the Cv value in function of  \n")+
                        string(" temperature T for the 2D Ising lattice. The test is performed \n")+
                        string(" for different lattice sizes: 2x2, 4x4, 8x8, 16x16. Critical \n")+
                        string(" temperature is Tc=2.27. From theoretical results we know that\n")+
                        string(" Cv should obtain maximum value near the Tc value. For large and\n")+
                        string(" very small values of T the Cv should tend to zero.\n")+
                        string(" The results are saved in following files:")+  
                        string(" 'IsingTestCC2D.txt' and 'IsingTestCC2D.png' in 'tests_out' directory. \n")+
                        string(" See run() function for more details.");
            
            info(); 
            run(); // run all calculations
            string cmd =string("cd ")+test_dir_output+string(";gnuplot IsingTestCC2D.plt");
            int info = system(cmd.c_str());
        }   
        
        /**
         * Tests Cv in function of T for different lattice sizes NxN for 2D Ising model
         */
        void run(){
         
        cout << " Running test..." << endl;                
        cout << " It may take some time..." << endl;
        
        //SIMULATION INPUT DATA
        int chainLength      =   2;    //initial lattice size
        double magneticField = 0.0;
        int therm_t          =  50;    //thermalization time
        int measure_f        =   1;    //measure frequency
        int prod_t           =1000;    //productiontime
        int maxGdist         =  10;    //domain of correlation function extends from 0 to maxGdist-1
        int maxTsep          =   2;
        double initState     = 0.7;    //values between 0.5 and 1.0
    
                
        vector<vector<double> > dataNxN[4];
        vector<double> dataT;
        // preparation of temperature array
        for(double T = 1.0 ; T < 5.0 ; T += 0.1 ){ 
            dataT.push_back(T);
        }
        // ---------------------------------------------------------
        // Calculation for NxN lattice.
        // ---------------------------------------------------------
        for(int i=0 ; i < 4 ; i++){                    
            chainLength  = pow(2,i+1); // choosing proper lattice size
            cout << "Starting simulation for lattice "<<chainLength << "x" <<  chainLength << "...\n";  
            
            for(unsigned int t = 0 ; t < dataT.size() ; t++ ){ 
                double T = dataT[t];
                cout << "T=" << T << "\t:";
                Ising2D lattice2d(chainLength, T, magneticField, initState, maxGdist, maxTsep);        
                lattice2d.MC_simulation(therm_t, prod_t, measure_f,METHOD_WOLFF);
                
                double error = lattice2d.ERROR(gettotalFname().c_str(), ERROR_CHI);
                double cc    = lattice2d.CC(gettotalFname().c_str());                
                vector<double> data;
                data.push_back(cc); data.push_back(error); dataNxN[i].push_back(data);
            }// end of for            
        } // end of for over lattice sizes
        
        // saving data to file
        string fileout = test_dir_output+"IsingTestCC2D.txt";
        ofstream data_out(fileout.c_str());     
        
        for(unsigned int t = 0 ; t < dataT.size() ; t++ ){
             data_out << std::scientific << dataT[t] << "\t";   
             for(int i = 0; i < 4 ; i++){
                 data_out << dataNxN[i][t][0] << "\t";
                 data_out << dataNxN[i][t][1] << "\t";
             }
             data_out << endl;
        }
                
        }// end of run()

};


#endif	/* ISINGTESTCC2D_H */

