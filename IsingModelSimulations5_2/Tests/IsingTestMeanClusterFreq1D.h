
#ifndef ISINGTESTMeanClusterFreq1D_H
#define	ISINGTESTMeanClusterFreq1D_H

#include "IsingTest.h"

class IsingTestMeanClusterFreq1D : public IsingTest {
public:
 
    IsingTestMeanClusterFreq1D() : IsingTest() {
        test_name = "Mean cluster size in function of temperature";
        test_info = string(" Run test to calculate the distribution of cluster size versus temperature \n") +
                string(" for the 1D Ising chains using Wolf MC algorithm") +
                string("\n The results are saved in proper files:") +
                string(" 'IsingTestMeanClusterFreq1D.txt' and 'IsingTestMeanClusterFreq1D.png' in 'tests_out' directory. \n") +
                string(" See run() function for more details.");

        info();
        run(); // run all calculations
        string cmd = string("cd ") + test_dir_output + string(";gnuplot IsingTestMeanClusterFreq1D.plt");
        int info = system(cmd.c_str());
    }

    /**
     * Tests cluster size in function of temperature for 1D Ising model
     */
    void run() {

        cout << " Running test..." << endl;

        //SIMULATION INPUT DATA
        int chainLength = 1000;
        double temperature = 2.0;
        double magneticField = 0.0;
        int therm_t = 10; //thermalization time
        int measure_f = 10; //measure frequency
        int prod_t = 1000; //productiontime
        int maxGdist = 10; //domain of correlation function extends from 0 to maxGdist-1
        int maxTsep = 2;
        double initState = 0.7; //values between 0.5 and 1.0
        string fileout = test_dir_output + "IsingTestMeanClusterFreq1D.txt";
        ofstream data_out(fileout.c_str());

        //CREATION AND INITIALIZATION OF THE INSTANCE OF ISING CLASS
        Ising chain(chainLength, temperature, magneticField, initState, maxGdist, maxTsep);
        Ising chain2(chainLength, temperature, magneticField, initState, maxGdist, maxTsep);

        //CALL TO THE ISING CLASS MEMBER FUNCTION PERFORMING SIMULATION
        //WOLF AND METROPOLIS CHAINS THAT ALLOWS TO COMPARE TWO ALGORITHMES
        chain.MC_simulation(therm_t, prod_t, measure_f,METHOD_WOLFF);
        chain2.MC_simulation(therm_t*10, prod_t*10, measure_f,METHOD_METROPOLIS);
        vector<vector<short> > SS;
        cluster_stat CS;
        float mean_cs_size,mean_cs_size2;
        for(int i=0 ; i<100 ; i++){
            if(i==99){
                chain.T=0;
                chain2.T=0;
            }
            //INITIALIZATION OF SPINS CHAINS
            SS.clear();
            for (int t = 0; t < prod_t; t++){
                if (!(t % measure_f)) {
                    SS.push_back(chain.S); 
                    chain.cycle();
                }else
                    chain.cycle();
            }
            
            //MEASUREMENT OF MEAN CLUSTER SIZE FOR GIVEN TEMPERATURE
            mean_cs_size=0;
            CS=chain.mean_cluster_freq1D(SS);
            for(unsigned int j=0; j< CS.CF.size(); j++){
                mean_cs_size+=CS.CF[j]*j;
            }
            //data_out << std::scientific << chain.T << "\t" << mean_cs_size << endl;
            chain.T-=0.02;
            
            //TERMALIZATION AFTER TEMPERATURE STEP
            for(int k=0;k<therm_t;k++) chain.cycle();
            
            //Part for metropolis chain
            //INITIALIZATION OF SPINS CHAINS
            SS.clear();
            for (int t = 0; t < prod_t*10; t++){
                if (!(t % measure_f*10)) {
                    SS.push_back(chain2.S); 
                    chain2.cycle();
                }else
                    chain2.cycle();
            }
            
            //MEASUREMENT OF MEAN CLUSTER SIZE FOR GIVEN TEMPERATURE
            mean_cs_size2=0;
            CS=chain2.mean_cluster_freq1D(SS);
            for(unsigned int j=0; j< CS.CF.size(); j++){
                mean_cs_size+=CS.CF[j]*j;
            }
            data_out << std::scientific << chain.T << "\t" << mean_cs_size << "\t" << mean_cs_size2 << endl;
            chain2.T-=0.02;
            
            //TERMALIZATION AFTER TEMPERATURE STEP
            for(int k=0;k<therm_t*10;k++) chain2.cycle();
        }
        data_out.close();
    }// end of run()

};

#endif	/* ISINGTESTMeanClusterFreq1D_H */

