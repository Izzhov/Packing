/*
 * NbrList2D.h
 *
 *  Created on: May 26, 2015
 *      Author: izzhov
 */

#ifndef NBRLIST2D_H_
#define NBRLIST2D_H_

#include "Torus2D.h"
#include "mymath.h"
#include <vector>
#include <cmath>

class NbrList2D {
private:
	Torus2D& torus;
	std::vector<std::vector<int> > nbrs;
	std::vector<std::vector<double> > locs;
	double L; //what L was when nbrlist was made
public:
	NbrList2D(Torus2D& torus);

	bool are_nbrs(int i, int j);
	bool should_rebuild();

	void build_nbrs();
	void build_locs();
	void build_L();
	void build_all();

	double get_L(){return L;}
};

#endif /* NBRLIST2D_H_ */
