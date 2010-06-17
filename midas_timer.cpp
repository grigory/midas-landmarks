// (c) Microsoft Corporation. All rights reserved.

/*************************************
 *
 * MIDASTimer (used for time measured)
 *
 *************************************/

#include "midas_timer.h"
#define INFINITE_TIME 10e100

/********************
 * Generic functions
 ********************/

MIDASTimer::MIDASTimer (bool s) {
	base_time = 0.0;
	max_time = 0.0;
	if (s) Start();
	else running = false;
}

MIDASTimer::MIDASTimer (double maxtime) {
	base_time = 0.0;
	SetMaxTime(maxtime);
	Start();
}

double MIDASTimer::GetTime() {
	if (running) return (GetElapsedTime() + base_time);
	else return base_time;
}

double MIDASTimer::Start() {
	double current = GetTime();
	base_time = 0.0;
	StartTiming();
	running = true;
	return current;
}

void MIDASTimer::SetMaxTime (double mt) {max_time = mt;}
double MIDASTimer::GetMaxTime () {return max_time;}

bool MIDASTimer::IsTimeExpired () {
	if (GetMaxTime() == 0) return false;
	bool time_expired = (GetTime() >= GetMaxTime());
	return time_expired;
}

double MIDASTimer::GetTimeToExpire () { //may be negative!
	if (GetMaxTime() == 0) return INFINITE_TIME;
	else return (GetMaxTime() - GetTime());
}

void MIDASTimer::SetBaseTime (double bt) {base_time = bt;}

double MIDASTimer::Reset() {
	double current = GetTime();
	running = false;
	base_time = 0.0;
	return current;
}

double MIDASTimer::Pause() {
	base_time = GetTime();
	running = false;
	return base_time;
}

double MIDASTimer::Resume() {
	if (running) return GetTime();
	else {
		running = true;
		StartTiming();
		return base_time;
	}
}


//-------------
// MIDAS_CLOCK
//-------------

void MIDASTimer::StartTiming() {start_time = GetUserTime();}

double MIDASTimer::GetElapsedTime() {
	return (GetUserTime() - start_time);
}

double MIDASTimer::GetUserTime() {
	double msecs = (double) clock() / CLOCKS_PER_SEC;
	if (msecs > 0) return msecs; 
	else return 0.0; //sometimes msecs is -0.000 (go figure...)
}
