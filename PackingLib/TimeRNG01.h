/*
 * TimeRNG01.h
 *
 *  Created on: May 22, 2015
 *      Author: root
 */

#ifndef TIMERNG01_H_
#define TIMERNG01_H_

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

class TimeRNG01 {
private:
	boost::random::mt19937 rng;
	boost::random::uniform_01<double> roll;
	int seed;
public:
	TimeRNG01();
	TimeRNG01(int t);//to seed directly

	double num();

	void set_seed(int seed);
	int get_seed(){return seed;}
};

#endif /* TIMERNG01_H_ */
