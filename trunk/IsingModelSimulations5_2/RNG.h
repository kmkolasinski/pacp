/* 
 * File:   RNG.h
 *
 */

#ifndef RNG_H
#define	RNG_H

//#########################################################################
//THIS FILE CONTAINS
//
//long int random_int(long int& j)
//
//double drandom(long int& seed2)
//
//void initializeRNG(long int& seed1)
//
//########################################################################



/************************************************************************
*Programs for Random Number Generation					*
* pages 451 - 452 of Newmann & Barkema book:				*
*"Monte Carlo Methods in Statistical Physics"	   		        *
*									*
* Program Descriptions: Three parts:					*
* A) Generator of random integers by linear congruency.  The generated	* 
*    integers are long int's in interval [0, 2^(31) - 2], i.e., this    *
*    gives 2 147 483 646 values.  See Newmann & Barkema, pages 386-389. *
* B) Shuffler of integers that produces random double floats in [0, 1).	*
* C) Initializer Program that must be called only once!!  This program	*
*	 seeds the Shuffler Program.					*	
*************************************************************************/

/* A) Generator of random integers by linear congruency.*/

// PRE-COND:	Calling program provides an initialized global, long int, seed 
//		parameter j.  Initial value of j must be in interval [1, m - 1].  
//		In program, m - 1 is selected to be: 2147483646 = 2^31 - 2.  
// POST-CONDITION:Returns positive random long integer in interval [0, m - 1].
//		Side effect is that seed parameter "j" is updated so that 
//		a sequence of calls to random_int(j) produces a random 
//		sequence of integers.  Here, "j" is updated, because "j" 
//		is a reference parameter in random_int(j).


//Declarations and initializations for random integer generator

class RNG {
public:
    



//GENERATOR LICZB PSEUDOLOSOWYCH FORWARD DECLARATIONS
static long int random_int(long int&);   // forward declaration -- random long int generator
static void initializeRNG(long int&);	 // forward declaration -- seeds random number generators
static double drandom(long int&);	 // forward declaration -- random double float generator


static const long int m; // m-1 = number of random integers generable
static const long int a; // values of constants a, q, and r give good 
static const long int q; // randomization and maximum cyle length m-1.
static const long int r;


//Additional declarations and initializations needed for shuffler program
static long int  jj;
static double denom;
static const int NR;
static long int    y;		// random int mixed into shuffling array.
static long int j[64];		// vector of random integers that are shuffled
			        // Both y and j[] are seeded by making the single call 
                           	// of initializeRNG(seed1).
};


#endif	/* RNG_H */

