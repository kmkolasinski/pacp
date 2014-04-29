/* 
 * File:   Tools.h
 *
 */

#ifndef TOOLS_H
#define	TOOLS_H


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

using namespace std;

void Display(vector<int> S) ;

vector<short> readDATAtoVector(string totalFname) ;

vector<vector<short> > readDATAtoVectors(string totalFname, int N);

void displayDATAfromVectors(string totalFname, int N);

string gettotalFname() ;


//ACCESS TO DATA STORED IN vector<short>
//    vector<short> simRes = readDATAtoVector1(totalFname);
//    vector<short>::iterator it;
//    for (it=simRes.begin(); it != simRes.end(); it++) cout<<" "<<*it;

//ACCESS TO DATA STORED IN vector<vector<short> > 





#endif	/* TOOLS_H */

