/*
 * Sizes.h
 *
 *  Created on: Jun 8, 2015
 *      Author: izzhov
 */

#ifndef SIZES_H_
#define SIZES_H_

#include <gsl/gsl_vector.h>

#include <vector>
#include <ctime>
#include <cmath>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>

class Sizes {
private:
	int N;
	std::vector<std::vector<gsl_vector*> > szs;
	std::vector<double> wts;
	double totwt;
	std::vector<std::vector<gsl_vector*> > fullsizes;
public:
	Sizes();
	Sizes(int N, std::vector<std::vector<gsl_vector*> > szs, std::vector<double> wts);

	//creates the vector fullsizes with sizes in correct proportions
	void populate();

	void set_N(int N);
	void set_szs(std::vector<std::vector<gsl_vector*> > szs);
	void set_wts(std::vector<double> wts);

	int get_N(){return N;}
	std::vector<std::vector<gsl_vector*> > get_szs(){return szs;}
	std::vector<double> get_wts(){return wts;}
	std::vector<std::vector<gsl_vector*> > get_fullsizes(){return fullsizes;}
};

#endif /* SIZES_H_ */
