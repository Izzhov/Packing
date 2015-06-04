/*
 * SpheroSizes.h
 *
 *  Created on: May 27, 2015
 *      Author: izzhov
 */

#ifndef SPHEROSIZES_H_
#define SPHEROSIZES_H_

#include <vector>
#include <ctime>
#include <cmath>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>

class SpheroSizes {
private:
	int N;
	std::vector<double> rs;
	std::vector<double> as;
	std::vector<double> wts;
	double totwt;
	std::vector<double> radii;
	std::vector<double> lengths; //actually alphas but w/e
public:
	SpheroSizes();
	SpheroSizes(int N, std::vector<double> rs, std::vector<double> as, std::vector<double> wts);

	void construct(int N, std::vector<double> rs, std::vector<double> as, std::vector<double> wts);

	void populate();//creates vector of radii and vector of alphas in correct prop's

	void set_N(int N);
	void set_rs(std::vector<double> rs);
	void set_as(std::vector<double> as);
	void set_wts(std::vector<double> wts);

	int get_N(){return N;}
	std::vector<double> get_rs(){return rs;}
	std::vector<double> get_as(){return as;}
	std::vector<double> get_wts(){return wts;}
	std::vector<double> get_radii(){return radii;}
	std::vector<double> get_lengths(){return lengths;}
};

#endif /* SPHEROSIZES_H_ */
