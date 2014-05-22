/* 
 * File:   IsingTestChi2D.h
 * Author: edyta
 *
 * Created on 22 maj 2014, 01:04
 */

#ifndef ISINGTESTCHI2D_H
#define	ISINGTESTCHI2D_H


#include "IsingTest.h"

class IsingTestChi2D : public IsingTest {

    public:
        IsingTestChi2D():IsingTest(){
            test_name = "Magnetic susceptibility test.";   
            test_info = string(" Run test to calculate the Chi value in function of \n")+
                        string(" temperature T for the 2D Ising matrix.\n The results are saved in proper files:")+  
                        string(" 'IsingTestChi2D.txt' and 'IsingTestChi2D.png' in 'tests_out' directory. \n")+
                        string(" See run() function for more details.\n");
                        string(" We plot the dependence of Chi on temperature. We expect that for low T");
                        string(" Chi explodes. We find a good agreement between the numerical and ");
                        string(" analytical solution (see IsingTestChi2D.png)");
            
            info(); 
            run(); // run all calculations
            string cmd =string("cd ")+test_dir_output+string(";gnuplot IsingTestChi2D.plt");
            int info = system(cmd.c_str());
        }   
        
        /**
         * Tests Chi in function of Temperature for 2D model
         */
        void run(){
         
        cout << " Running test..." << endl;    
            
        //SIMULATION INPUT DATA
        int chainLength      = 5;        
        double magneticField = 0.0;
        int therm_t          = 100;   //thermalization time
        int measure_f        = 1;     //measure frequency
        int prod_t           = 1000;  //productiontime
        int maxGdist         = 10;     //domain of correlation function extends from 0 to maxGdist-1
        int maxTsep          = 2;
        double initState     = 0.7;    //values between 0.5 and 1.0
        string fileout = test_dir_output+"IsingTestChi2D.txt";
        string fileout2 = test_dir_output+"IsingTestChi2D_2.txt";
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
        Ising2D chain2D(chainLength, temperature, magneticField, initState, maxGdist, maxTsep);
        //CALL TO THE ISING CLASS MEMBER FUNCTION PERFORMING SIMULATION
        chain2D.MC_simulation(therm_t, prod_t, measure_f, METHOD_METROPOLIS);
                
        double chi      = chain2D.Chi(gettotalFname().c_str());
        //double chi_anal = exp(2./temperature)/temperature;
        double error    = chain2D.ERROR(gettotalFname().c_str(), ERROR_CHI);

        // saving data to file
        data_out << std::scientific <<  temperature  << "\t" 
                                         //<< chi_anal << "\t"
                                         << chi      << "\t"
                                         << error    << "\t"
                                         << chi-error << "\t"
                                         << chi+error << endl;

        cout << std::scientific <<  temperature  << "\t" << chi << "\t" /*<< chi_anal */<< "\t" <<  error << endl;

        }// end of for(temp)
        data_out.close();     
        
        
        ofstream data_out2(fileout2.c_str());
        double temperature=1.;
        for( int n=2 ; n < 21 ; n += 1  ){         
        //CREATION AND INITIALIZATION OF THE INSTANCE OF ISING CLASS
        Ising2D chain2D(n, temperature, magneticField, initState, maxGdist, maxTsep);
        //CALL TO THE ISING CLASS MEMBER FUNCTION PERFORMING SIMULATION
        chain2D.MC_simulation(therm_t, prod_t, measure_f, METHOD_METROPOLIS);
                
        double chi      = chain2D.Chi(gettotalFname().c_str());
        //double chi_anal = exp(2./temperature)/temperature;
        double error    = chain2D.ERROR(gettotalFname().c_str(), ERROR_CHI);

        // saving data to file
        data_out2 << std::scientific <<  n  << "\t" 
                                         //<< chi_anal << "\t"
                                         << chi      << "\t"
                                         << error    << "\t"
                                         << chi-error << "\t"
                                         << chi+error << endl;

        cout << std::scientific <<  temperature  << "\t" << chi << "\t" /*<< chi_anal */<< "\t" <<  error << endl;

        }// end of for(temp)
        data_out2.close(); 
        
        
        
        }// end of run()

};


#endif	/* ISINGTESTCHI2D_H */

