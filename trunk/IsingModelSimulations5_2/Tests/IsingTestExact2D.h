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
            test_name = "Exact solution for Ising model: exact specific heat";   
            test_info = string(" Run test to calculate the exact value of specific heat Cv for 2D lattice. \n")+
                        string(" The results are saved in proper files:\n")+  
                        string(" 'IsingTestExact2D.txt' and 'IsingTestExact2D.png' in 'tests_out' directory. \n")+
                        string(" In docs folder you can find: 1) a chapter of book ('IsingTestExactBook.pdf') in which you can find the explanation of the \n")+
                        string(" method we used to computed exat value of Cv, and 2) the short theoretical explanation is provided in 'IsingExactSolution.pdf' file. \n")+
                        string(" In brief, we run loop over all possible \n")+
                        string(" configurations for given N by N lattice using the Gray code enumeration method. We calculate\n")+
                        string(" the density of states N(E) - number of confs. with the same energy E. Then we use \n")+
                        string(" N(E) to calculate partition function Z = Sum_E N(E) exp(-beta*E), from which we calculate\n")+
                        string(" Cv from definition using standard numerical differentiation.\n")+
                        string(" This method allows us to get exact solution for 2D lattice. But is limited to very small lattices,\n")+
                        string(" up to 5x5.\n")+
                        string(" You may compare result obtained from this method with IsingTestCC2D.png file.\n")+
                        string(" which was obtained from MC simulation.\n")+
                        string(" See run() function for more details.");
            
            info(); 
            run(); // run all calculations
            string cmd =string("cd ")+test_dir_output+string(";gnuplot IsingTestExact2D.plt");
            int info = system(cmd.c_str());
        }   
        
        void run(){
         
        cout << " Running test..." << endl;    
        // Here we run calculation from 2 to 5. Where nx means the size of the lattice Nx by Nx.
        for(nx = 2; nx<6; nx++){
        cout << "Calculating exact heat capacity for lattice: " << nx << "x" << nx << endl;    
        N  = nx*nx;
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
            NE[E()] += 1;
            S[k] *= -1;
        }    
        // normalization of N(E) - dividing by number of states
        for(it = NE.begin() ; it != NE.end() ; it++){
            it->second = it->second / pow(2,N);                                
        } 
        
        stringstream ss;
        ss << nx;
        string fileout = test_dir_output+"IsingTestExact2D"+ss.str()+"x"+ss.str()+".txt";
        ofstream data_out(fileout.c_str());
        for(T = 0.1 ; T < 10 ; T += 0.1){ // temperature dependence loop
            double dT = 0.01;
            double C1 = log(Z(T+dT,NE));
            double C2 = log(Z(T-dT,NE));
            double C3 = log(Z(T,NE));
            double Cv = 2*T*( C1 - C2 )/2/dT + T*T*( C1 + C2 - 2*C3 )/dT/dT;
            // here we calculate Cv using numerical derivative, but we choose small
            // value of dT which should gives us still accurate results
            data_out << T << "\t" << Cv/N << endl; // Cv per spin
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
};

#endif	/* ISINGTESTEXACT2D_H */

