#include "RNG.h"

const long int RNG::m=2147483647; // m-1 = number of random integers generable
const long int RNG::a=16807;		 // values of constants a, q, and r give good 
const long int RNG::q=127773;	 // randomization and maximum cyle length m-1.
const long int RNG::r=2836;


//Additional declarations and initializations needed for shuffler program
long int  RNG::jj =89343;
double    RNG::denom=(1.0/(m-1));
const int RNG::NR = 64;

long int  RNG::y = 0;
long int  RNG::j[64] = {0};	// vector of random integers that are shuffled
			        // Both y and j[] are seeded by making the single call 
                           	// of initializeRNG(seed1).
//_____________________________________________________________________________________________

long int RNG::random_int(long int& j)
{
	long int l;
	l = j/q;	// Parameter "j" is recursively defined to produce a
			// a sequence of random numbers.  
	j = a*(j - q*l) -r*l;	// This iteration step is actually a mod 2^31 
				// operation as described by Schrage.  Here, first term = 
				// j mod q and second term is (int divide(j/q) )*r.  
				// See Newmann & Barkema, page 388.
	if (j<0) j += m;	// This shift bring the number back into the interval
							// [1, m - 1].
	return (j-1);	// j - 1 is a random integer in interval [0, m - 2]
}	                //end program to return random long integer




/* B) Remainder for Shuffler program that generates double float in [0,1). */

// PRE-COND:	Variable seed2 has been declared and initialized in [1, 2^31 - 2]. 
//		Also, a single call must be made to initializeRNG(seed1) to seed the 
//		y and j[] variables of the Shuffler program.  

// POST-COND:	Shuffler returns random double float in [0, 1).  Repeat cycle is 
//		believed to be about (m - 1)^64 where m = 2^31 - 1.  Thus, m is 
//		about = 2x10^9, and the repeat cycle, i.e., (m - 1)^64 is thus, 
//		more or less infinite in any reasonably doable simulation.
//		Including this file in a program makes each call of 
//		drandom(long int& i) update i, yy, and j[k] as a side effect.  The
//		reason is: seed2, y, and j[] are global variables.
 	

//_____________________________________________________________________________________________

double RNG::drandom(long int& seed2)
{
	long int l;
	long int k;

	l = seed2/q;
	seed2 = a*(seed2 - q*l) - r*l;	//Selects "seed2" to be a random integer in 
			// [1, 2^31 - 2 ] and updates "seed2" as reference variable.
			//This ensures that a series of calls will produce a
			//a series of random long ints "seed2".
	if (seed2 < 0)	//ensures that seed2 is random int in [1, m - 1]
		{
			seed2 = seed2 + m;
		}
	k = (y/m)*NR;		//Selects an integer in [0, NR]

	y = j[k];		// Sets y = integer of k-th component of j[].
				// Side effect: global variable y is updated for 
				// next call of drandom.
	j[k] = seed2;		// Side effect: global variable j[k] is updated 
				// for next call of drandom.

	return denom*(y - 1);	// Returns randomly shuffled double float in [0, 1).
}	//end Shuffler program



/* C) Program to initialize Shuffler program. */

// PRE-COND:	Variable seed1 has been declared and initialized in [1, 2^31 - 2].  
//				Note: variable y and vector j[] are globally declared in a program 
//				that includes this header file. 
// POST-COND:	Shuffler variable y and components of vector j[] are seeded with 
//				randomly selected long integers in [1, m - 1].
//				Side effect of fact that seed1 is reference variable is that 
//				each initialization of a component of j[] updates seed1 thereby 
//				ensuring that the sequence of components are given a random sequence  
//				of long int values.

//_____________________________________________________________________________________________

void RNG::initializeRNG(long int& seed1) {
    int s;
    for (s = 0; s < NR; s++) //loop iteratively seeds each component of j[].
    {
        j[s] = random_int(seed1); //seeds component "s" of shuffler vector j[].
        //	std::cout << "j[" << s << "] = " << j[s] << '\n'; 
        //	The above line is used as a test to insure that j[] is a random array
        //	That is, to print out array j[] to check random initialization	
    }
    y = random_int(seed1); // seeds variable j
    //	std::cout << "yy = " << yy<< '\n'; 
    //	That is, to print out yy to check random initialization
} //end program to initialize Shuffler program



