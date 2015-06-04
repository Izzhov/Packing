/*
 * TimeRNG01.cpp
 *
 *  Created on: May 22, 2015
 *      Author: root
 */

#include "TimeRNG01.h"

TimeRNG01::TimeRNG01() {
	unsigned int t = static_cast<unsigned int>(time(0));
	rng.seed(t);
	set_seed(t);
}

TimeRNG01::TimeRNG01(int t){
	rng.seed(t);
	set_seed(t);
}

double TimeRNG01::num() {
	return roll(rng);
}

void TimeRNG01::set_seed(int seed){
	this->seed = seed;
}
