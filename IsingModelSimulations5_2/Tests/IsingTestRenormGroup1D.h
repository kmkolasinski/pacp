/* 
 * File:   IsingTestRenormGroup1D.h
 * Author: mkk
 *
 * Created on 5 de Maio de 2014, 18:09
 */

#ifndef ISINGTESTRENORMGROUP1D_H
#define	ISINGTESTRENORMGROUP1D_H


#include "IsingTest.h"

class IsingTestRenormGroup1D : public IsingTest {

   public:
       IsingTestRenormGroup1D():IsingTest(){
           test_name = "Test of renormalisation group transformation with majority rule.";   
           test_info = string(" Run test to get dependence of temperature (calculated from \n")+
                       string(" correlation) on step of renormalisation.\n")+
                       string(" We start with lattice of length which is power of 3.\n")+
                       string(" We use majority rule to assign a spin for each block of 3 spins\n")+
                       string(" and calculate temperature (derived from correlations) for such a shorter chain.\n")+
                       string(" Then we repeat few times the whole operation.\n")+   
                       string(" We expect that for shorter chains the temperature will wander off the critical point:\n")+
                       string(" in 1-dimensional case it will go further form T=0 K.\n")+
                       string(" \n The results are saved in proper files:")+  
                       string(" 'IsingTestRenormGroup1D.txt' and 'IsingTestRenormGroup.png' in 'tests_out' directory. \n")+
                       string(" See run() function for more details.");
           
           
           
           info(); 
           run(); // run all calculations
           string cmd =string("cd ")+test_dir_output+string(";gnuplot IsingTestRenormGroup1D.plt");
           int info = system(cmd.c_str());
       }   
       
       /**
        * Tests Chi in function of Temperature for 1D model
        */
       void run(){
        
         cout << " Running test..." << endl;    
         
         int chainLength      = 100;        
         double magneticField = 0.0;
         int therm_t          = 100;   //thermalization time
         int measure_f        = 1;     //measure frequency
         int prod_t           = 10000;  //productiontime
         int maxGdist         = 20;     //domain of correlation function extends from 0 to maxGdist-1
         int maxTsep          = 2;
         double initState     = 0.7;    //values between 0.5 and 1.0    
         double temp          = 1.0;
         
         string fileout = test_dir_output+"IsingTestRenormGroup1D.txt";
         ofstream data_out(fileout.c_str());
         
         int coeff=81;
        
         cout << "\n \n chainLength = " << chainLength*coeff << endl;
         
        Ising chain3(chainLength*coeff,temp,magneticField,initState,maxGdist,maxTsep);   
        chain3.MC_simulation(therm_t,prod_t,measure_f);
        chain3.mean_S_correlation(gettotalFname());
        
                data_out << std::scientific <<  chainLength*coeff  << "\t" 
                                 << chain3.Gs[1] << "\t" 
                                 << (1.0 / atanh(chain3.Gs[1] / chain3.Gs[0])) << endl;

      
         
        for(coeff=81;coeff>1;){ 
        
        vector<vector<short> > SS3 = readDATAtoVectors(gettotalFname(), chainLength*coeff);

        vector<vector<short> >::iterator t;
        vector<short>::iterator pos;
        
        fstream DATA(gettotalFname().c_str(), ios::out);
        DATA.setf(ios::showpos);
        
        for (unsigned int i =  0; i < SS3.size() ; i++) { 
            for (int k =  0; k < SS3[0].size() ; k+=3 ) { //loop over chain sites "pos" at fixed time
                DATA <<  (( SS3[i][k] + SS3[i][k+1] + SS3[i][k+2]  ) >0?1:-1) << " ";
                //cout << i << ":" << SS3[i][k]  <<  "  " <<  SS3[i][k+1] << "  " <<  SS3[i][k+2] ;
                //cout << (( SS3[i][k] + SS3[i][k+1] + SS3[i][k+2]  ) >0?1:-1) << endl;
            } 
//            for (unsigned int k =  0; k < SS3[0].size() ; k+=3 ) { 
//                DATA << SS3[i][k] << " ";
//                DATA << SS3[i][k] << " ";
//                DATA << SS3[i][k] << " ";
//            } 
            DATA << endl;
        }
        DATA.close();
        coeff/=3;
        cout << "\n \n chainLength = " << chainLength*coeff << endl;
        Ising chain1(chainLength*coeff,temp,magneticField,initState,maxGdist,maxTsep);   
        chain1.mean_S_correlation(gettotalFname()); 
        
                data_out << std::scientific <<  chainLength*coeff  << "\t" 
                                         << chain1.Gs[1] << "\t" 
                                         << (1.0 / atanh(chain1.Gs[1] / chain1.Gs[0])) << endl;

        
        }
        
        data_out.close();
                                      
       }// end of run()

};


#endif	/* ISINGTESTRENORMGROUP1D_H */

