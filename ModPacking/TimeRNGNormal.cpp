/*
 * TimeRNGNormal.cpp
 *
 *  Created on: Jun 22, 2015
 *      Author: izzhov
 */

#include "TimeRNGNormal.h"

TimeRNGNormal::TimeRNGNormal() {
	unsigned int t = static_cast<unsigned int>(time(0));
	rng.seed(t);
	set_seed(t);
}

TimeRNGNormal::TimeRNGNormal(int t) {
	rng.seed(t);
	set_seed(t);
}

double TimeRNGNormal::num() {
	return roll(rng);
}

void TimeRNGNormal::set_seed(int seed) {
	this->seed = seed;
}
