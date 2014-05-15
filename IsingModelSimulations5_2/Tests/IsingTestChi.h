
#ifndef ISINGTESTCHI_H
#define	ISINGTESTCHI_H

#include "IsingTest.h"

class IsingTestChi : public IsingTest {

    public:
        IsingTestChi():IsingTest(){
            test_name = "Magnetic susceptibility test.";   
            test_info = string(" Run test to calculate the Chi value in function of \n")+
                        string(" temperature T for the 1D Ising chain.\n The results are saved in proper files:")+  
                        string(" 'IsingTestChi1D.txt' and 'IsingTestChi1D.png' in 'tests_out' directory. \n")+
                        string(" See run() function for more details.\n");
                        string(" We plot the dependence of Chi on temperature. We expect that for low T");
                        string(" Chi explodes. We find a good agreement between the numerical and ");
                        string(" analytical solution (see IsingTestChi.png)");
            
            info(); 
            run(); // run all calculations
            string cmd =string("cd ")+test_dir_output+string(";gnuplot IsingTestChi1D.plt");
            int info = system(cmd.c_str());
        }   
        
        /**
         * Tests Chi in function of Temperature for 1D model
         */
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
        string fileout = test_dir_output+"IsingTestChi1D.txt";
        ofstream data_out(fileout.c_str());
        
        // the range of the temperature : (0.3,2) 
        cout << " Start  T=" << 0.3 << endl;
        cout << " End    T=" << 2.0 << endl;
        cout << " Step  dT=" << 0.1 << endl;
        cout << " It may take some time..." << endl;
        cout << " T           " << " Chi(numerical) " << " Chi(analytical)" << "    Std. dev" << endl;
        cout << "----------------------------------------------------" << endl;
        for( double temperature = 0.3 ; temperature < 2.0 ; temperature += 0.1  ){         
        //CREATION AND INITIALIZATION OF THE INSTANCE OF ISING CLASS
        Ising chain(chainLength, temperature, magneticField, initState, maxGdist, maxTsep);
        //CALL TO THE ISING CLASS MEMBER FUNCTION PERFORMING SIMULATION
        chain.MC_simulation(therm_t, prod_t, measure_f);
                
        double chi      = chain.Chi(gettotalFname().c_str());
        double chi_anal = exp(2./temperature)/temperature;
        double error    = chain.ERROR(gettotalFname().c_str(), ERROR_CHI);

        // saving data to file
        data_out << std::scientific <<  temperature  << "\t" 
                                         << chi_anal << "\t"
                                         << chi      << "\t"
                                         << error    << "\t"
                                         << chi-error << "\t"
                                         << chi+error << endl;

        cout << std::scientific <<  temperature  << "\t" << chi << "\t" << chi_anal << "\t" <<  error << endl;

        }// end of for(temp)
        data_out.close();                                    
        }// end of run()

};


#endif	/* ISINGTESTCHI_H */

