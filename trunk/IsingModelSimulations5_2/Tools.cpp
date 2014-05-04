/* 
 * File:   Tools.cpp
 */


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

using namespace std;

void Display(vector<int> S) {
    for (int i = 0; i < S.size(); i++) {
        (S[i] == 1) ? cout << " o" : cout << " x";
    }
}

vector<short> readDATAtoVector(string totalFname) {
    //Domain: tools, Ising model in 1D.
    //Read simulation data from disk file and put them to vector.
    //If N represents the length of the spin chain then each bloc of N elements  
    //of the vector represents a snapshot of the state of system.
    //totalFname is the string "path/filename"
    fstream DATA(totalFname.c_str(), ios::in);
    vector<short> S;
    int s;
    while (DATA >> s)
        S.push_back(s);
    return S;
}

vector<vector<short> > readDATAtoVectors(string totalFname, int N) {
    //Domain: tools, Ising model in 1D.
    //Read simulation data from disk file and put them to vector S of vectors s
    //Vector s of size = length of the spin chain represents a state of the system.
    //The size of the vector S is equal to number of measures.
    //The elements of S are states of system ordered chronologically
    fstream DATA(totalFname.c_str(), ios::in);
    int i = 0, j = 0, x;
    vector<vector<short> > SS;
    vector<short> s;
    while (DATA >> x) {
        s.push_back(x);
        i++;
        if (!(i % N)) {
            SS.push_back(s);
            s.clear();
        };
    }
    return SS;
}

void displayDATAfromVectors(string totalFname, int N) {
    vector<vector<short> > SS = readDATAtoVectors(totalFname, N);
    vector<vector<short> >::iterator t; //time iterator
    vector<short> ss;
    vector<short>::iterator pos;        //space iterator

    for (t = SS.begin(); t != SS.end(); t++) {
        for (pos = t->begin(); pos != t->end(); pos++) cout << " " << *pos;
        //Display states of the chain in lines; each lines contain instantaneous 
        //state of the chain; lines are ordered chronologically.
        cout << endl;
    }
}

string gettotalFname() {
    string path, fname, totalFname; //path to the project folder on your laptop
    //***Replace this path by the path on your laptop
    path = "";
    fname = "MCdata";
    totalFname = path + fname;
    return totalFname;
}



