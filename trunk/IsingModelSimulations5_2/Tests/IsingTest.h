#ifndef ISINGTEST_H
#define	ISINGTEST_H

#include <iostream>
#include <cstdlib>
#include <string>
#include "../Ising2D.h"
#include "../Tools.h"



using namespace std;
class IsingTest {//this is abstract class

public:
    
    IsingTest(){
        test_dir_output = "tests_out/";        
    };
    
    /**
     *  Print info about simulation.
     */ 
    void info(){        
        cout << "----------------------------------------------------" << endl;
        cout << " Starting test: " << test_name << endl;
        cout << "----------------------------------------------------" << endl;
        cout << test_info << endl;        
        cout << "----------------------------------------------------" << endl;
    };
    
    /**
     * Start abstract test need to be implemented.
     */ 
    virtual void run() = 0;
    
protected:
    string test_name;       // name of test
    string test_dir_output; // output directory
    string test_info;       // info message
};

#endif	/* ISINGTEST_H */

