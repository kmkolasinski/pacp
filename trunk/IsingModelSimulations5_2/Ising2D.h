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
#define Yfrom2D(site)( int(site/nx) ) // map from 1D array to y coordinate of 2D array
#define Xfrom2D(site)(  (site%nx) )   // map from 1D array to x coordinate of 2D array

struct cluster_stat_2D {
    vector<float> CF; //CF[x]=nb of occurences of the spin 
                      //cluster of length x in a given state of the Ising chain 
    
    vector<short> SS; //Additional vector allowing to create image of not only spins
                      //but also clusters.
    int cmax; //maximal length of the cluster in the distribution
} ;

class Ising2D : public Ising {
public:
    Ising2D();
    Ising2D(int NN,double TT, double hh, double p, int maxGd, int maxTsep); //Constructor of the class
    virtual ~Ising2D();   
    double E(); 
    cluster_stat_2D cluster_freq2D();
//*** MONTE CARLO ENGINS, MC WRAPPER, SOME MC PARAMETERS     
    vector<short> Metropolis_cycle();
    vector<short> Wolff_cycle();
    void MC_simulation(int therm_t, int prod_t, int measure_f,ISING_METHOD_TYPE mt = METHOD_METROPOLIS);



    int nx; //size of box
};

#endif	/* ISING2D_H */

