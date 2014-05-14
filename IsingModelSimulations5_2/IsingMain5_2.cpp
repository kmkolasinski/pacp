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

using namespace std;

int main(int argc, char** argv) {
    
    // ----------------------------------------------------------
    // Comment or uncomment following tests:
    // ----------------------------------------------------------

    // IsingTestChi isingTestChi;     //  Ising1D test of Chi in function of T
    // IsingTestError isingTestError; //  Ising1D test of ERROR funtion for Chi for different prod_
    // IsingTestClusterStat isingTestClusterStat;  //Ising1D test of cluster distribution
    // IsingTestTCorr isingTestTCorr; // 
    // IsingTestCC2D isingTestCC2D;   //  Ising2D test of CC in function of T for different lattices
    // IsingTestM2D isingTestM2D;     //  Ising2D test of M in function of T for different lattices        
    // IsingTestExact2D isingTestExact2D; // Ising2D test which calculates the exacts values of Cv for  small lattices
    // IsingTestRenormGroup1D  isingTestRenormGroup1D; // Ising1D test for renormalization group

    
    cout << "\n\n     T H E   E N D " << endl;
    return 0;
}

