/* 
* File:   IsingTestError.h
* Author: alina
*
* Created on 4 maj 2014, 18:11
*/

#ifndef ISINGTESTERROR1D_H
#define	ISINGTESTERROR1D_H

#include "IsingTest.h"

class IsingTestError1D : public IsingTest {

   public:
       IsingTestError1D():IsingTest(){
           test_name = "Test of error (bootstrap method).";   
           test_info = string(" Run test to calculate the error value in function of \n")+
                       string(" the production time using the bootstrap method. For more details see:\n ")+  
                       string(" http://en.wikipedia.org/wiki/Bootstrapping_(statistics).\n ")+  
                       string(" \nThe results are saved in following files:")+  
                       string(" 'IsingTestError.txt' and 'IsingTestError1D.png' in 'tests_out' directory. \n")+
                       string(" This test shows that the error decreases with number of prod_t. We founded that \n")+
                       string(" the asymptotic behavior of error(prod_t) behaves like 1/sqrt(prod_t) which \n")+
                       string(" agrees with the general result for error known in MC integration methods. For more \n")+
                       string(" details see: http://en.wikipedia.org/wiki/Monte_Carlo_integration. \n")+
                       string(" See run() function for more details.");
           
           info(); 
           run(); // run all calculations
           string cmd =string("cd ")+test_dir_output+string(";gnuplot IsingTestError1D.plt");
           int info = system(cmd.c_str());
       }   
       
       /**
        * Tests Error for 1D chain
        */
       void run(){
        
       cout << " Running test..." << endl;    
           
       //SIMULATION INPUT DATA
       int chainLength      = 500;        
       double magneticField = 0.0;
       int therm_t          = 500;   //thermalization time
       int measure_f        = 1;     //measure frequency

       int maxGdist         = 10;     //domain of correlation function extends from 0 to maxGdist-1
       int maxTsep          = 2;
       double initState     = 0.7;    //values between 0.5 and 1.0
       string fileout = test_dir_output+"IsingTestError1D.txt";
       string fileout2 = test_dir_output+"IsingTestError1D_T.txt";
       ofstream data_out(fileout.c_str());
       
       cout << " Start prod_t=" << 2000    << endl;
       cout << " End   prod_t=" << 500000  << endl;
       cout << " Mult  prod_t=" << 1.2       << endl;
       cout << " It may take some time..." << endl;
       cout << " prod_t   " << " Error "   << endl;
       cout << "----------------------------------------------------" << endl;
       
       double temperature = 1.0;
       Ising chain(chainLength, temperature, magneticField, initState, maxGdist, maxTsep);
       for ( int prod_t = 2000 ; prod_t < 50000 ; prod_t *= 1.2  ){
           chain.MC_simulation(therm_t, prod_t, measure_f);
           double error=chain.ERROR(gettotalFname().c_str(), ERROR_CHI);
           double chi=chain.Chi(gettotalFname().c_str());
           data_out << std::scientific <<  prod_t << "\t" << error <<  "\t" << chi << "\n";
           cout     << std::scientific <<  prod_t << "\t" << error <<  "\t" << chi << "\n";
       }// end of for(prod_t)
       
       data_out.close();  
       
       
       ofstream data_out2(fileout2.c_str());
       int prod_t = 5000;
       for ( double temperature2 = 0.5 ; temperature2 < 3.05 ; temperature2 += 0.1  ){
           Ising chain(chainLength, temperature2, magneticField, initState, maxGdist, maxTsep);
           chain.MC_simulation(therm_t, prod_t, measure_f);
           double error=chain.ERROR(gettotalFname().c_str(), ERROR_CHI);
           double chi=chain.Chi(gettotalFname().c_str());
           double chi_a=exp(2./temperature2)/temperature2/chain.N;
           data_out2 << std::scientific <<  temperature2 << "\t" << error <<  "\t" << chi <<  "\t" << chi_a << "\n";
           cout     << std::scientific <<  temperature2 << "\t" << error <<  "\t" << chi <<  "\t" << chi_a << "\n";
       }// end of for(prod_t)
       
       data_out2.close();  
       
       
       
       
       
       
       
       
       
       
       
       }// end of run()

};


#endif	/* ISINGTESTERROR1D_H */
