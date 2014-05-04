/* 
 * File:   Ising2D.h
 * Author: mkk
 *
 * Created on 4 de Maio de 2014, 9:13
 */

#ifndef ISING2D_H
#define	ISING2D_H

#include "Ising5_2.h"

#define To2D(i,j)( i + (j)*nx )

class Ising2D : public Ising {
public:
    Ising2D();
    Ising2D(int NN,double TT, double hh, double p, int maxGd, int maxTsep); //Constructor of the class
    virtual ~Ising2D();   
    double E(); 
    
//*** MONTE CARLO ENGINS, MC WRAPPER, SOME MC PARAMETERS     
    vector<short> Metropolis_cycle();
    void MC_simulation(int therm_t, int prod_t, int measure_f);


private:
    int nx; //size of box
};

#endif	/* ISING2D_H */

