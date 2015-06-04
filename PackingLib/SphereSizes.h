/*
 * SphereSizes.h
 *
 *  Created on: May 21, 2015
 *      Author: root
 */

#ifndef SPHERESIZES_H_
#define SPHERESIZES_H_

#include <vector>
#include <ctime>
#include <cmath>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>

class SphereSizes {
private:
	int N;
	std::vector<double> sizes;
	std::vector<double> wts;
	double totwt;
	std::vector<double> radii;
public:
	SphereSizes();
	SphereSizes(int N, std::vector<double> sizes, std::vector<double> wts);

	void cons(int N, std::vector<double> sizes, std::vector<double> wts);

	void construct(int N, std::vector<double> sizes, std::vector<double> wts);

	void populate();

	void set_N(int N);
	void set_sizes(std::vector<double> sizes);
	void set_wts(std::vector<double> wts);

	int get_N(){return N;}
	std::vector<double> get_sizes(){return sizes;}
	std::vector<double> get_wts(){return wts;}
	std::vector<double> get_radii(){return radii;}
};

#endif /* SPHERESIZES_H_ */
