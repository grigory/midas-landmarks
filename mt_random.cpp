// (c) Microsoft Corporation. All rights reserved.

#include "mt_random.h"

const unsigned long int MTRandom::maxvalue = (unsigned long int)(-1); 

const unsigned long int MTRandom::MATRIX_A = 0x9908b0dfUL;   /* constant vector a */
const unsigned long int MTRandom::UPPER_MASK = 0x80000000UL; /* most significant w-r bits */
const unsigned long int MTRandom::LOWER_MASK = 0x7fffffffUL; /* least significant r bits */

unsigned long MTRandom::mt[N]; /* the array for the state vector  */
int MTRandom::mti=N+1; /* mti==N+1 means mt[N] is not initialized */
			
