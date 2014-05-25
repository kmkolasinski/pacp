
#ifndef ISINGTESTCORRELATIONLENGTH_H
#define	ISINGTESTCORRELATIONLENGTH_H

#include "IsingTest.h"

class IsingTestCorrelationLength : public IsingTest {
public:
    
    //Additional informations about fit appear while running the program due to use of 
    //function print in gnuplot. Each correlation length is obtained from its plot versus
    //correlation distance.
    
    IsingTestCorrelationLength() : IsingTest() {
        test_name = "Comparison between mean cluster size and correlation length";
        test_info = string(" Run test to calculate the distribution of domains of parallel spins \n") +
                string(" for the 1D Ising chain and for fixed temperature using Wolf MC algorithm") +
                string("\n The results are saved in proper files:") +
                string(" 'IsingTestCLvsMCS.txt' and 'IsingTestCLvsMCS.png' in 'tests_out' directory. \n") +
                string(" See run() function for more details.");

        info();
        run(); // run all calculations
        string cmd = string("cd ") + test_dir_output + string(";gnuplot IsingTestCorrelationLength.plt");
        int info = system(cmd.c_str());
        //ADDITIONAL OPERATIONS ON PLOT AND DATA FILES
        cmd = string("cd ") + test_dir_output + string(";paste -d ' ' cl.txt IsingTestCorrelationLengthMCS.txt > IsingTestCLvsMCS.txt");
        info=system(cmd.c_str());
        cmd = string("cd ") + test_dir_output + string(";gnuplot IsingTestCorrelationLengthB.plt;clear");
        info=system(cmd.c_str());
        //;rm cl.txt fit.log
    }

    /**
     * Tests relation between correlation length and mean cluster size.
     */
    void run() {

        cout << " Running test..." << endl;

        //SIMULATION INPUT DATA
        int chainLength = 1000;
        double temperature = 1.5;
        double magneticField = 0.0;
        int therm_t = 100; //thermalization time
        int measure_f = 10; //measure frequency
        int prod_t = 1000; //productiontime
        int maxGdist = 10; //domain of correlation function extends from 0 to maxGdist-1
        int maxTsep = 2;
        double initState = 0.7; //values between 0.5 and 1.0
        string fileout = test_dir_output + "IsingTestCorrelationLengthMCS.txt";
        ofstream data_out(fileout.c_str());
        fileout = test_dir_output + "IsingTestCorrelationLength.txt";
        ofstream data_out_cl(fileout.c_str());
        
        //CREATION AND INITIALIZATION OF THE INSTANCE OF ISING CLASS
        Ising chain(chainLength, temperature, magneticField, initState, maxGdist, maxTsep);
        
        //FIRST MC SIMULATION USING WOLF ALGORITHM
        chain.MC_simulation(therm_t, prod_t, measure_f,METHOD_WOLFF);
        
        //TEMPORARY STRUCTURES CONTAINING SPIN CHAINS IN SS AND MEAN SPIN CORRELATION.
        //ALL NEEDED CALCULATIONS OF CORRELATION LENGTHS ARE HANDLEN BY GNUPLOT
        //THANKS TO IT'S USEFUL FITTING OPTIONS.
        
        vector<vector<short> > SS;
        cluster_stat CS;
        vector <vector<double> > GS;
        float mean_cs_size;
        for(int i=0;i<10;i++){
            //CREATION OF SPIN CHAINS IN GIVEN TEMPERATURE
            SS.clear();
            for (int j=0; j<prod_t; j++){
                if (!(j % measure_f)) {
                    SS.push_back(chain.S);

                    chain.cycle();
                }else
                    chain.cycle();
            }   
            mean_cs_size=0;
            CS=chain.mean_cluster_freq1D(SS);
            for(unsigned int j=0; j< CS.CF.size(); j++){
                mean_cs_size+=CS.CF[j]*j;
            }
            data_out << std::scientific << mean_cs_size << endl;
            
            //ADDING NEW CORRELATION VECTOR
            GS.push_back(chain.mean_S_correlation(SS));
            
            chain.T-=0.1;
            for(int k=0;k<therm_t;k++) chain.cycle();
        }
        data_out.close();
        for(unsigned int i=0 ; i<GS[0].size() ; i++){
            data_out_cl << i;
            for(unsigned int j=0 ; j<GS.size() ; j++){
                data_out_cl << "\t" << GS[j][i];
            }
            data_out_cl << endl;
        }
        data_out_cl.close();
    }// end of run()

};

#endif	/* ISINGTESTCORRELATIONLENGTH_H */

