/*
 * TimeRNGNormal.h
 *
 *  Created on: Jun 22, 2015
 *      Author: izzhov
 */

#ifndef TIMERNGNORMAL_H_
#define TIMERNGNORMAL_H_

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

class TimeRNGNormal {
private:
	boost::random::mt19937 rng;
	boost::random::normal_distribution<double> roll;
	int seed;
public:
	TimeRNGNormal();
	TimeRNGNormal(int t);//to seed directly

	double num();

	void set_seed(int seed);
	int get_seed(){return seed;}
};

#endif /* TIMERNGNORMAL_H_ */
