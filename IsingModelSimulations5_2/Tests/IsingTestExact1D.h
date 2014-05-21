/* 
 * File:   IsingTestExact2D.h
 * Author: edyta
 *
 * Created on 8 maj 2014, 20:01
 */

#ifndef ISINGTESTEXACT1D_H
#define	ISINGTESTEXACT1D_H



#include "IsingTest.h"
#include <map>

class IsingTestExact1D : public IsingTest , public Ising {

    public:
        IsingTestExact1D():IsingTest(){
            test_name = "Exact solution for Ising model in 1D: exact specific heat";   
            test_info = string(" Run test to calculate the exact value of specific heat Cv for 1D chain. \n")+
                        string(" The results are saved in proper files:")+  
                        string(" 'IsingTestExact1DN=x.txt' (where x denotes chain lenght) and 'IsingTestExact1D.png' in 'tests_out' directory. \n")+
                        string(" We also put a chapter of book ('IsingTestExact2D.pdf') in which you can find the explanation of the \n")+
                        string(" method we used to computed exat value of Cv. Short explanation of code is located in 'IsingTestExactComments.pdf' file. \n")+
                        string(" In brief, we run loop over all possible \n")+
                        string(" configurations for given chain length N using the Gray code enumeration method. We calculate\n")+
                        string(" the density of states N(E) - number of confs. with the same energy E. Then we use \n")+
                        string(" density of states - N(E) to calculate partition function Z = Sum_E N(E) exp(-beta*E), from which we calculate\n")+
                        string(" Cv from definition using standard numerical differentiation.\n")+
                        string(" This method allows us to get exact solution for 1D chain. But is limited to very small chains up to N=25.\n")+    
                        string(" See run() function for more details.");
            
            info(); 
            run(); // run all calculations
            string cmd =string("cd ")+test_dir_output+string(";gnuplot IsingTestExact1D.plt");
            int info = system(cmd.c_str());
        }   
        
        void run(){
         
        cout << " Running test..." << endl;    
        // Here we run calculation from N=2 to 22 (2,7,12,17,22). Where N means the chain length.
        for( N = 2 ; N < 26 ; N+=5 ){
        cout << "Calculating exact heat capacity for chain N=" << N << endl;    
        
        h  =   0.0; 
        S.clear();
        // We create initial state
        for (int i = 0; i < N; i++) {
            S.push_back(-1);
        }
        
        // t - is auxilliary vector used by gray flip procedure
        vector<int> t(N+1);
        for( int i = 1 ; i < N+2 ; i++){
            t[i-1] = i;
        }
        map<int,double> NE; // density of states N(E)
        map<int,double>::iterator it;
        
        // calculation of N(E) - loop over all possible states
        for(int i=0; i < pow(2,N);i++){
            int k = gray_flip(t);
            int iE = E()*N;
            NE[iE]  +=  1;            
            S[k]    *= -1;
        }    
        // normalization of N(E) - dividing by number of states        
        for(it = NE.begin() ; it != NE.end() ; it++){            
            it->second = it->second / pow(2,N);                                
        } 
        
        stringstream ss;
        ss << N;
        string fileout = test_dir_output+"IsingTestExact1DN="+ss.str()+".txt";
        ofstream data_out(fileout.c_str());
        for(T = 0.1 ; T < 10 ; T += 0.01){ // temperature dependence loop
            double dT = 0.01;
            double C1 = log(Z(T+dT,NE));
            double C2 = log(Z(T-dT,NE));
            double C3 = log(Z(T,NE));
            double Cv = 2*T*( C1 - C2 )/2/dT + T*T*( C1 + C2 - 2*C3 )/dT/dT;
            // here we calculate Cv using numerical derivative, but we choose small
            // value of dT which should gives us still accurate results
            data_out << T << "\t" << Cv/N << "\t" << exactCv(T,N) << endl; // Cv per spin
        }        
        data_out.close();        
        }// end of nx loop              
        }// end of run()
        
        
        // -----------------------------------------------------------------
        // Auxilliary functions
        // -----------------------------------------------------------------
        /**
         * Performs gray flip on current lattice state S. Compare this with
         * Algorithm 5.2 from tests_out/IsingTestExact2D.pdf
         */ 
        int gray_flip(vector<int>& t){
            int k = t[0];
            if( k > N ) return k-1;
            t[k-1] = t[k];
            t[k] = k+1;
            if(k!=1) t[0] = 1;
            return k-1;
        };
        
        /**
         * Caluculates the partition function from definition. Compare this 
         * with section 5.1.2 of tests_out/IsingTestExact2D.pdf 
         */ 
        double Z(double T,map<int,double>& NE){
            double Z = 0;
            map<int,double>::iterator it;
            for(it = NE.begin() ; it != NE.end() ; it++){
                Z += it->second * exp( -it->first / T );               
            }          
            return Z;
        }
        /**
         * Analytical value of specific heat Cv for 1D chain and B = 0 T. 
         * @param T - temperature
         * @param N - chain length
         * @return specific heat Cv of 1D chain
         */
        double exactCv(double T,int N){
            
            double dT = 0.001;
            // Here we perform numerical derivation:
            double Z  = 4*pow(cosh(1/T),N)*(pow(tanh(1/T),N)+1);
            double Zp = 4*pow(cosh(1/(T+dT)),N)*(pow(tanh(1/(T+dT)),N)+1);
            double Zm = 4*pow(cosh(1/(T-dT)),N)*(pow(tanh(1/(T-dT)),N)+1);
            
            return T*( (T+dT)*log(Zp) + (T-dT)*log(Zm) - 2*(T)*log(Z) )/dT/dT/N;
        }

};








#endif	/* ISINGTESTEXACT1D_H */

