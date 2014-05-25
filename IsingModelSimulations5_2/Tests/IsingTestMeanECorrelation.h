
#ifndef ISINGTESTMEANECORRELATION_H
#define	ISINGTESTMEANECORRELATION_H

#include "IsingTest.h"

class IsingTestMeanECorrelation : public IsingTest {
public:

    IsingTestMeanECorrelation() : IsingTest() {
        test_name = "Mean energy correlation test function";
        test_info = string(" Run test to calculate the mean energy correlation distribution versus corr distance \n") +
                string(" for the 1D Ising chains using Wolf MC algorithm") +
                string("\n The results are saved in proper files:") +
                string(" 'IsingTestMeanECorrelation.txt' and 'IsingTestMeanECorrelation.png' in 'tests_out' directory. \n") +
                string(" See run() function for more details.");

        info();
        run(); // run all calculations
        string cmd = string("cd ") + test_dir_output + string(";gnuplot IsingTestMeanECorrelation.plt");
        int info = system(cmd.c_str());
    }

    /**
     * Tests distribution of mean energy correlation versus correlation
     * distance in given temperature.
     */
    void run() {

        cout << " Running test..." << endl;

        //SIMULATION INPUT DATA
        int chainLength = 1000;
        double temperature = 1.0;
        double magneticField = 0.0;
        int therm_t = 100; //thermalization time
        int measure_f = 100; //measure frequency
        int prod_t = 10000; //productiontime
        int maxGdist = 20; //domain of correlation function extends from 0 to maxGdist-1
        int maxTsep = 2;
        double initState = 0.7; //values between 0.5 and 1.0
        string fileout = test_dir_output + "IsingTestMeanECorrelation.txt";
        ofstream data_out(fileout.c_str());

        //CREATION AND INITIALIZATION OF THE INSTANCE OF ISING CLASS
        Ising chain(chainLength, temperature, magneticField, initState, maxGdist, maxTsep);
        //CALL TO THE ISING CLASS MEMBER FUNCTION PERFORMING SIMULATION
        chain.MC_simulation(therm_t, prod_t, measure_f);
    
        vector<vector<short> > SS;
        //INITIALIZATION OF SPINS CHAINS
        for(int j=0 ; j<10 ; j++ ){    
            SS.clear();
            for (int t=0; t<prod_t ; t++ ){
                if (!(t % measure_f)) SS.push_back(chain.S); 
                chain.cycle();
            }
        }
        cout << "Spin chains initialization complete." << endl;
        vector<double> mEc;
        mEc=chain.mean_E_correlation(SS);
        
        for(int i=0 ; i<mEc.size() ; i++){
            data_out<< i << "\t" << mEc[i] << endl; 
        }
        data_out.close();
    }// end of run()

};

#endif	/* ISINGTESTMeanClusterFreq1D_H */

