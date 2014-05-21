/* 
 * File:   IsingMain5_2.cpp
 *
 */

//////////////////////////////////////////////////////////////////////////////
// MAIN FUNCTION
//////////////////////////////////////////////////////////////////////////////

#include <stdio.h>      
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>		
#include <cstring>		
#include <iomanip>
#include <cmath>

//#include<stack> //LIFO (Last In First Out) container for Wolff MC algorithm

#include "Tools.h"
#include "Ising5_2.h"
#include "Ising2D.h"
#include "Tests/IsingTestChi.h"
#include "Tests/IsingTestClusterStat.h"
#include "Tests/IsingTestCC2D.h"
#include "Tests/IsingTestError.h"
#include "Tests/IsingTestM2D.h"
#include "Tests/IsingTestExact2D.h"
#include "Tests/IsingTestRenormGroup1D.h"
#include "Tests/IsingTestTCorr.h"
#include "Tests/IsingTestMeanClusterFreq1D.h"
#include "Tests/IsingTestCorrelationLength.h"
#include "Tests/IsingTestMeanECorrelation.h"
#include "Tests/IsingTestE1D.h"

using namespace std;

int main(int argc, char** argv) {

//CALCULATION OF PHYSICAL PROPERTIES OF ISING MODEL IN ONE AND TWO DIMENSIONS 
//                 IN FUNCTION OF THE TEMPERATURE. 
//       Periodic boundary conditions are assumed in all cases.    
//Special attention is paid to behavior close to criticality i.e. to T=0 in 1D
//and to T=2.27 for 2D square lattice.
//To run any of the following test please uncomment it, compile and run the program.
//To modify input data edit test and change the value of parameters, then compile and run.    
    
    // ----------------------------------------------------------
    // Comment or uncomment following tests:
    // ----------------------------------------------------------
   
//1D ISING CHAIN WITH PERIODIC BOUNDARY CONDITIONS 
    // IsingTestChi isingTestChi;     //  Ising1D test of Chi in function of T
    // IsingTestError isingTestError; //  Ising1D test of ERROR funtion for Chi for different prod_
    // IsingTestClusterStat isingTestClusterStat;  //Ising1D test of cluster distribution
    // IsingTestCorrelationLength isingTestCorrelationLength;
     IsingTestMeanECorrelation isingTestMeanECorrelation;
    // IsingTestTCorr isingTestTCorr; // 
    // IsingTestMeanClusterFreq1D isingTestMeanClusterFreq1D; // Ising1D test for mean cluster size distribution
    // IsingTestRenormGroup1D  isingTestRenormGroup1D; // Ising1D test for renormalization group
    
//2D ISING MODEL ON SQUARE LATTICE WITH PERIODIC BOUNDARY CONDITIONS 
    // IsingTestCC2D isingTestCC2D;   //  Ising2D test of CC in function of T for different lattices
    // IsingTestM2D isingTestM2D;     //  Ising2D test of M in function of T for different lattices        
    // IsingTestExact2D isingTestExact2D; // Ising2D test which calculates the exacts values of Cv for  small lattices
    // IsingTestE1D isingTestE1D;

//TOOLS
    // TestError TestError; //  Generic test of ERROR funtion for some input data 
    //                           and any algorithm serving for calculation of physical
    //                           quantity.
        
    cout << "\n\n     T H E   E N D " << endl;
    return 0;
}