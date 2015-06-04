/*
 * SpheroTorus2D.h
 *
 *  Created on: May 27, 2015
 *      Author: izzhov
 */

#ifndef SPHEROTORUS2D_H_
#define SPHEROTORUS2D_H_

#include "Sphero2D.h"
#include "SpheroSizes.h"
#include "TimeRNG01.h"
#include "mymath.h"
#include <vector>
#include <cmath>
#include <algorithm>

class SpheroTorus2D {
private:
	int N;
	double L;
	SpheroSizes ssizes;
	TimeRNG01 trng;
	std::vector<Sphero2D> sphs;
	double partarea;
public:
	SpheroTorus2D();
	SpheroTorus2D(int N, std::vector<double> rs,
			std::vector<double> as, std::vector<double> wts);
	SpheroTorus2D(int N, double L, std::vector<double> rs,
			std::vector<double> as, std::vector<double> wts);
	//makes alpha=2, r=1 parts in set locs
	SpheroTorus2D(int N, double L, std::vector<std::vector<double> > locs);

	//populates randomly
	void populate();
	//fills with pre-set locations
	void set_populate(std::vector<std::vector<double> > locs);

	void set_N(int N);
	void set_L(double L);
	void set_1x(int i, double x);
	void set_1y(int i, double y);
	void set_u(int i, double u);
	void set_v(int i, double v);

	int get_N(){ return N; }
	double get_L(){ return L; }

	double get_1x(int i){return sphs.at(i).get_x();}
	double get_1y(int i){return sphs.at(i).get_y();}
	double get_u(int i){return sphs.at(i).get_u();}
	double get_v(int i){return sphs.at(i).get_v();}

	double get_x(int i){return L*mymath::pmod(sphs.at(i).get_x(),1);}
	double get_y(int i){return L*mymath::pmod(sphs.at(i).get_y(),1);}
	double get_r(int i){return sphs.at(i).get_r();}
	double get_a(int i){return sphs.at(i).get_a();}

	void find_set_L(double phi);

	void normalize(int i);

	double rel_x(int i, int j, int k); //range [0,1)
	double rel_y(int i, int j, int k);
	double len(int i); //length of bone (line segment)
	double near(int i, int j); //smallest ctr-ctr dist that ever matters
	double R(int i, int j); //smallest bone-bone dist that ever matters
	double Iom(int i);

	double rij2(int i, int j, int k);//algorithm: see papers-to-cite
	double ui_rij(int i, int j, int k);
	double uj_rij(int i, int j, int k);
	double ui_uj(int i, int j, int k);

	double lip(int i, int j, int k); //algorithm: see papers-to-cite
	double ljp(int i, int j, int k);
	double bli(int i, int j, int k);
	double blj(int i, int j, int k);
	double lis(int i, int j, int k);
	double ljs(int i, int j, int k);
	double ell2(int i, int j, int k); //shortest dist. btwn bones, squared

	//packing fraction
	double pack_frac();
};

#endif /* SPHEROTORUS2D_H_ */
