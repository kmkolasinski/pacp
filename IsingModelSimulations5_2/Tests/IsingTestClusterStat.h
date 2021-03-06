
#ifndef ISINGTESTCLUSTERSTAT_H
#define	ISINGTESTCLUSTERSTAT_H

#include "IsingTest.h"

class IsingTestClusterStat : public IsingTest {
public:

    IsingTestClusterStat() : IsingTest() {
        test_name = "Distribution of spin domains";
        test_info = string(" Run test to calculate the distribution of domains of parallel spins \n") +
                string(" for the 1D Ising chain and for fixed temperature using Metropolis MC algorithm") +
                string("\n The results are saved in proper files:") +
                string(" 'IsingTestClusterStat.txt' and 'IsingTestClusterStat.png' in 'tests_out' directory. \n") +
                string(" See run() function for more details.");

        info();
        run(); // run all calculations
        string cmd = string("cd ") + test_dir_output + string(";gnuplot IsingTestClusterStat.plt");
        int info = system(cmd.c_str());
    }

    /**
     * Tests domains distribution in function of temperature for 1D model
     */
    void run() {

        cout << " Running test..." << endl;

        //SIMULATION INPUT DATA
        int chainLength = 1000;
        double temperature = 1.5;
        double magneticField = 0.0;
        int therm_t = 1000; //thermalization time
        int measure_f = 10; //measure frequency
        int prod_t = 10000; //productiontime
        int maxGdist = 10; //domain of correlation function extends from 0 to maxGdist-1
        int maxTsep = 2;
        double initState = 0.7; //values between 0.5 and 1.0
        string fileout = test_dir_output + "IsingTestClusterStat.txt";
        ofstream data_out(fileout.c_str());

        //CREATION AND INITIALIZATION OF THE INSTANCE OF ISING CLASS
        Ising chain(chainLength, temperature, magneticField, initState, maxGdist, maxTsep);
        Ising chain2(chainLength, temperature, magneticField, initState, maxGdist, maxTsep);
        //CALL TO THE ISING CLASS MEMBER FUNCTION PERFORMING SIMULATION
        //WOLF AND METROPOLIS CHAINS THAT ALLOWS TO COMPARE TWO ALGORITHMES
        chain.MC_simulation(therm_t/10, prod_t/100, measure_f,METHOD_WOLFF);
        chain2.MC_simulation(therm_t, prod_t, measure_f);
        
        //cluster_stat clStat = chain.cluster_freq1D(gettotalFname().c_str());
        cluster_stat clStat = chain.cluster_freq();
        cluster_stat clStatb = chain2.cluster_freq();
        for (int L = 1; L <= clStat.cmax; L++) {
            data_out << std::scientific << L << "\t" << clStat.CF[L] << "\t" << clStatb.CF[L] << endl;
        }

        data_out.close();
    }// end of run()

};

#endif	/* ISINGTESTCLUSTERSTAT_H */

