
#ifndef ISINGTESTE1D_H
#define	ISINGTESTE1D_H

#include "IsingTest.h"

class IsingTestE1D : public IsingTest {
public:

    IsingTestE1D() : IsingTest() {
        test_name = "Energy temperature distribution.";
        test_info = string(" Run test to calculate the energy distribution versus temperature \n") +
                string(" for the 1D Ising chains using Wolf MC algorithm") +
                string("\n The results are saved in proper files:") +
                string(" 'IsingTestE1D.txt' and 'IsingTestE1D.png' in 'tests_out' directory. \n") +
                string(" See run() function for more details.");

        info();
        run(); // run all calculations
        string cmd = string("cd ") + test_dir_output + string(";gnuplot IsingTestE1D.plt");
        int info = system(cmd.c_str());
    }

    /**
     * Tests distribution of mean energy correlation versus correlation
     * distance in given temperature.
     */
    void run() {

        cout << " Running test..." << endl;

        //SIMULATION INPUT DATA
        int chainLength = 500;
        double temperature = 0.0;
        double magneticField = 0.0;
        int therm_t = 10; //thermalization time
        int measure_f = 10; //measure frequency
        int prod_t = 1000; //productiontime
        int maxGdist = 20; //domain of correlation function extends from 0 to maxGdist-1
        int maxTsep = 2;
        double initState = 0.9; //values between 0.5 and 1.0
        string fileout = test_dir_output + "IsingTestE1D.txt";
        ofstream data_out(fileout.c_str());

        //CREATION AND INITIALIZATION OF THE INSTANCE OF ISING CLASS
        Ising chain(chainLength, temperature, magneticField, initState, maxGdist, maxTsep);
        //CALL TO THE ISING CLASS MEMBER FUNCTION PERFORMING SIMULATION
        chain.MC_simulation(therm_t, prod_t, measure_f,METHOD_WOLFF);

        //INITIALIZATION OF SPINS CHAINS
        float energy=0;
        for(int i=0 ; i<150 ; i++){    
            energy=chain.E();
            data_out << chain.T << "\t" << energy << endl;
            chain.T+=0.03;
            for(int j=0;j<therm_t;j++) chain.cycle();
        }
        data_out.close();
    }// end of run()

};

#endif	/* ISINGTESTE1D_H */

