// (c) Microsoft Corporation. All rights reserved.

/**********************************************************
 *
 * MIDASTimer: measures running time in seconds
 *
 **********************************************************/

#ifndef MIDAS_TIMER_H
#define MIDAS_TIMER_H

//-----------------------------------------------------
// heuristic to determine which timing facility to use
//-----------------------------------------------------

#include <stdio.h>
#include <time.h>

//------------------
// the class itself
//------------------
class MIDASTimer {
	private:
		bool running; //is it running now? (false -> paused)
		double base_time; //time of previous runs since last reset
		double max_time;  //reference time  

		
		double GetUserTime();
		double start_time;

		void SetBaseTime (double bt); 

		//facility-dependent functions
		void StartTiming();      //store time for future comparison
		double GetTime(); //return current time

		double Pause();   //pause and return current time
		double Resume();  //continue if paused, start if reset; return time before resuming
		double Start();   //reset and resume; return time before reset
		double Reset();   //reset timer and pause (at zero); return time before reset

		//auxiliary functions
		void SetMaxTime (double mt);

	public:
		//basic functions
		MIDASTimer (double maxtime);
		MIDASTimer (bool start=false);
		double GetElapsedTime(); //time since last resume/start
		double GetMaxTime();
		double GetTimeToExpire();

		bool IsTimeExpired ();
};

#endif
