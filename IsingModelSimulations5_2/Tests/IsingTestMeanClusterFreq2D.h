
#ifndef ISINGTESTMeanClusterFreq2D_H
#define	ISINGTESTMeanClusterFreq2D_H

#include "IsingTest.h"

class IsingTestMeanClusterFreq2D : public IsingTest {
public:
 
    IsingTestMeanClusterFreq2D() : IsingTest() {
        test_name = "Mean cluster size in function of temperature";
        test_info = string(" Run test to calculate the distribution of cluster size versus temperature \n") +
                string(" for the 2D Ising chains using Wolf MC algorithm") +
                string("\n The results are saved in proper files:") +
                string(" 'IsingTestMeanClusterFreq2D.txt' and 'IsingTestMeanClusterFreq2D.png' in 'tests_out' directory. \n") +
                string(" See run() function for more details.");

        info();
        run(); // run all calculations
        string cmd = string("cd ") + test_dir_output + string(";gnuplot IsingTestMeanClusterFreq2D.plt");
        int info = system(cmd.c_str());
    }

    /**
     * Tests cluster size in function of temperature for 1D Ising model
     */
    void run() {

        cout << " Running test..." << endl;

        //SIMULATION INPUT DATA
        int chainLength      = 20;    //initial lattice size
        double magneticField = 0.0;
        double T             = 5.0;
        int therm_t          =  50;    //thermalization time
        int measure_f        =   1;    //measure frequency
        int prod_t           =1000;    //productiontime
        int maxGdist         =  10;    //domain of correlation function extends from 0 to maxGdist-1
        int maxTsep          =   2;
        double initState     = 0.7;    //values between 0.5 and 1.0
 
        string fileout = test_dir_output+"IsingTestMeanClusterFreq2D.txt";
        ofstream data_out(fileout.c_str()); 
                
        Ising2D lattice_w(chainLength, T, magneticField, initState, maxGdist, maxTsep);
        //Ising2D lattice_m(chainLength, T, magneticField, initState, maxGdist, maxTsep);
        
        //WOLF AND METROPOLIS CHAINS THAT ALLOWS TO COMPARE TWO ALGORITHMES
        lattice_w.MC_simulation(therm_t, prod_t, measure_f,METHOD_WOLFF);
        //lattice_m.MC_simulation(therm_t, prod_t, measure_f);
       
        vector<vector<short> > SS;
        cluster_stat_2D CS;
        float mean_cs_size;
        for(int i=0 ; i<100 ; i++){
            cout << lattice_w.T << endl;
            //INITIALIZATION OF SPINS CHAINS
            SS.clear();
            for (int t = 0; t < prod_t; t++){
                if (!(t % measure_f)) {
                    SS.push_back(lattice_w.S); 
                    lattice_w.cycle();
                }else
                    lattice_w.cycle();
            }
            
            //MEASUREMENT OF MEAN CLUSTER SIZE FOR GIVEN TEMPERATURE
            mean_cs_size=0;
            CS=lattice_w.mean_cluster_freq2D(SS);
            for(unsigned int j=0; j< CS.CF.size(); j++){
                mean_cs_size+=CS.CF[j]*j;
            }
            data_out << std::scientific << lattice_w.T << "\t" << mean_cs_size << endl;
            lattice_w.T-=0.05;
            
            //TERMALIZATION AFTER TEMPERATURE STEP
            for(int k=0;k<therm_t;k++) lattice_w.cycle();
        }
        data_out.close();
    }// end of run()

};

#endif	/* ISINGTESTMeanClusterFreq1D_H */

