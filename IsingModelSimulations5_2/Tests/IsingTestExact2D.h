/* 
 * File:   IsingTestExact2D.h
 * Author: edyta
 *
 * Created on 8 maj 2014, 20:01
 */

#ifndef ISINGTESTEXACT2D_H
#define	ISINGTESTEXACT2D_H



#include "IsingTest.h"
#include <map>

class IsingTestExact2D : public IsingTest , public Ising2D {

    public:
        IsingTestExact2D():IsingTest(){
            test_name = "Exact solution for Ising model.";   
            test_info = string(" Run test to calculate the Chi value in function of \n")+
                        string(" temperature T for the 1D Ising chain.\n The results are saved in proper files:")+  
                        string(" 'IsingTestExact2D.txt' and 'IsingTestExact2D.png' in 'tests_out' directory. \n")+
                        string(" See run() function for more details.");
            
            info(); 
            run(); // run all calculations
            string cmd =string("cd ")+test_dir_output+string(";gnuplot IsingTestExact2D.plt");
            int info = system(cmd.c_str());
        }   
        
        /**
         * Tests Chi in function of Temperature for 1D model
         */
        void run(){
         
        cout << " Running test..." << endl;    
        for(nx = 2; nx<6; nx++){
            
        N  = nx*nx;
        h = 0;
        for (int i = 0; i < N; i++) {
            S.push_back(-1);
        }
        
        
        vector<int> t(N+1);
        for( int i = 1 ; i < N+2 ; i++){
            t[i-1] = i;
        }
        map<int,double> NE; // gestosc stanow
        map<int,double>::iterator it;
        
        for(int i=0; i < pow(2,N);i++){
            int k = gray_flip(t);
            NE[E()] += 1;
            S[k] *= -1;
            //if( (i%(pow(2,N)/10)) == 0) cout << float(i)/pow(2,N) << endl;
        }    
        int z_sum = 0;
        for(it = NE.begin() ; it != NE.end() ; it++){
            //cout << it->first << "  = " << it->second << endl;            
            z_sum += it->second;            
        } 
        for(it = NE.begin() ; it != NE.end() ; it++){
            it->second = it->second / z_sum;
            //cout << it->first << "  = " << it->second << endl;            
             
        } 
        
        stringstream ss;
        ss << nx;
        string fileout = test_dir_output+"IsingTestExact2D"+ss.str()+"x"+ss.str()+".txt";
        ofstream data_out(fileout.c_str());
        for(T = 0.1 ; T < 10 ; T += 0.1){
            double dT = 0.01;
            double C1 = log(partition_function(T+dT,NE));
            double C2 = log(partition_function(T-dT,NE));
            double C3 = log(partition_function(T,NE));
            double Cv = 2*T*( C1 - C2 )/2/dT + T*T*( C1 + C2 - 2*C3 )/dT/dT;
            data_out << T << "\t" << Cv/N << endl;
        }
        
        data_out.close();
        
        }
        
        
        }// end of run()
        
        int gray_flip(vector<int>& t){
            int k = t[0];
            if( k > N ) return k-1;
            t[k-1] = t[k];
            t[k] = k+1;
            if(k!=1) t[0] = 1;
            return k-1;
        };
        
        
        double partition_function(double T,map<int,double>& NE){
            double Z = 0;
            map<int,double>::iterator it;
            for(it = NE.begin() ; it != NE.end() ; it++){
                Z += it->second * exp( -it->first / T );               
            }          
            return Z;
        }

};








#endif	/* ISINGTESTEXACT2D_H */

