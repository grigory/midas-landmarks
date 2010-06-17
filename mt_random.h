// (c) Microsoft Corporation. All rights reserved.

/* C++ Wrapper written by Renato Werneck on original C code by 
   Nishimura and Matsumoto. Original header follows. */

/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

#ifndef MT_RANDOM_H
#define MT_RANDOM_H

#include <stdlib.h>
#include <stdio.h>

class MTRandom {
	private:
		static const unsigned long int maxvalue;

		/* Period parameters */  
		enum {N=624, M=397};
		static const unsigned long int MATRIX_A;   /* constant vector a */
		static const unsigned long int UPPER_MASK; /* most significant w-r bits */
		static const unsigned long int LOWER_MASK; /* least significant r bits */

		static unsigned long mt[N]; /* the array for the state vector  */
		static int mti; /* mti==N+1 means mt[N] is not initialized */

			
		/* initializes mt[N] with a seed */
		static void init_genrand(unsigned long s)
		{
			mt[0]= s & 0xffffffffUL;
			for (mti=1; mti<N; mti++) {
				mt[mti] = 
				(1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
				/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
				/* In the previous versions, MSBs of the seed affect   */
				/* only MSBs of the array mt[].                        */
				/* 2002/01/09 modified by Makoto Matsumoto             */
				mt[mti] &= 0xffffffffUL;
				/* for >32 bit machines */
			}
		}


		/* generates a random number on [0,0xffffffff]-interval */
		static unsigned long genrand_int32(void)
		{
			unsigned long y;
			static unsigned long mag01[2]={0x0UL, MATRIX_A};
			/* mag01[x] = x * MATRIX_A  for x=0,1 */

			if (mti >= N) { /* generate N words at one time */
				int kk;

				if (mti == N+1)   /* if init_genrand() has not been called, */
					init_genrand(5489UL); /* a default initial seed is used */

				for (kk=0;kk<N-M;kk++) {
					y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
					mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
				}
				for (;kk<N-1;kk++) {
					y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
					mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
				}
				y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
				mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

				mti = 0;
			}
		  
			y = mt[mti++];

			/* Tempering */
			y ^= (y >> 11);
			y ^= (y << 7) & 0x9d2c5680UL;
			y ^= (y << 15) & 0xefc60000UL;
			y ^= (y >> 18);

			return y;
		}

	public:
		//constructors
		MTRandom () {Randomize(1);}

		//randomize procedures
		static void Randomize (unsigned long s) {
			if (s==0) s = 1;
			init_genrand(s);
		}
		
		static unsigned long GetRand() {return genrand_int32();}

		//pick an integer uniformly at random between inf and sup (both inclusive)
		static int GetInteger (int inf, int sup) {
			if (sup<=inf) return inf;
			unsigned long range, minallowed, u;

			range = (unsigned long)(sup-inf+1);    //number of values allowed
			minallowed = (maxvalue % range) + 1;   //restrict search space to avoid small numbers
			if (minallowed==range) minallowed = 0;
			do {u = GetRand();}     //repeat until a good number is found
			while (u < minallowed);

			return (inf + (int)(u % range)); //return a number in the range
		}



		static float GetFloat () {return (float)GetDouble();} //get a float number in [0;1]
		static double GetDouble() {return GetDoubleOpen();} //double in the range [0;1]
		static double GetDoubleClosed() {return ((double)GetRand()/(double)maxvalue);}  //double in the range [0;1]
		static double GetDoubleOpen() {return ((double)GetRand()/((double)(maxvalue)+1.0));} //double in the range [0;1)
		static bool GetBool () {return (GetRand() & 1);}
};

#endif
