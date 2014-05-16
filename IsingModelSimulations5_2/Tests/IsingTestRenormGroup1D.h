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
           test_name = "Test of zly opis.";   
           test_info = string(" Run test to calculate the erro value in function of \n")+
                       string(" the production time.\n The results are saved in proper files:")+  
                       string(" 'IsingTestError.txt' and 'IsingTestError.png' in 'tests_out' directory. \n")+
                       string(" See run() function for more details.");
           
           
           
           info(); 
           run(); // run all calculations
          // string cmd =string("cd ")+test_dir_output+string(";gnuplot IsingTestError.plt");
          // int info = system(cmd.c_str());
       }   
       
       /**
        * Tests Chi in function of Temperature for 1D model
        */
       void run(){
        
         cout << " Running test..." << endl;    
         
         int chainLength      = 500;        
         double magneticField = 0.0;
         int therm_t          = 100;   //thermalization time
         int measure_f        = 1;     //measure frequency
         int prod_t           = 10000;  //productiontime
         int maxGdist         = 20;     //domain of correlation function extends from 0 to maxGdist-1
         int maxTsep          = 2;
         double initState     = 0.7;    //values between 0.5 and 1.0    
         double temp          = 1.4;
         
        Ising chain3(chainLength*3,temp,magneticField,initState,maxGdist,maxTsep);   
        chain3.MC_simulation(therm_t,prod_t,measure_f);
        chain3.mean_S_correlation(gettotalFname());        
//      
        vector<vector<short> > SS3 = readDATAtoVectors(gettotalFname(), chainLength*3);

        vector<vector<short> >::iterator t;
        vector<short>::iterator pos;
        int i = 0 ;
        fstream DATA(gettotalFname().c_str(), ios::out);
        DATA.setf(ios::showpos);
        
        for (int i =  0; i < SS3.size() ; i++) { 
//            for (int k =  0; k < SS3[0].size() ; k+=3 ) { //loop over chain sites "pos" at fixed time
//                DATA <<  (( SS3[i][k] + SS3[i][k+1] + SS3[i][k+2]  ) >0?1:-1) << " ";
//                //cout << i << ":" << SS3[i][k]  <<  "  " <<  SS3[i][k+1] << "  " <<  SS3[i][k+2] ;
//                //cout << (( SS3[i][k] + SS3[i][k+1] + SS3[i][k+2]  ) >0?1:-1) << endl;
//            } 
            for (int k =  0; k < SS3[0].size() ; k+=3 ) { 
                DATA << SS3[i][k] << " ";
//                DATA << SS3[i][k] << " ";
//                DATA << SS3[i][k] << " ";
            } 
            DATA << endl;
        }
        DATA.close();
        Ising chain1(chainLength,temp,magneticField,initState,maxGdist,maxTsep);   
        chain1.mean_S_correlation(gettotalFname());            
                              
       }// end of run()

};


#endif	/* ISINGTESTRENORMGROUP1D_H */

