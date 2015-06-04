/*
 * Torus2D.h
 *
 *  Created on: May 21, 2015
 *      Author: root
 */

#ifndef TORUS2D_H_
#define TORUS2D_H_

#include "Sphere2D.h"
#include "SphereSizes.h"
#include "TimeRNG01.h"
#include "mymath.h"
#include <vector>
#include <cmath>

class Torus2D {
private:
	int N;
	double L;
	SphereSizes ssizes;
	TimeRNG01 trng;
	std::vector<Sphere2D> sphs;
	double partarea;
	int seed;
public:
	Torus2D();
	Torus2D(int N, std::vector<double> sizes, std::vector<double> wts);
	Torus2D(int N, std::vector<double> sizes, std::vector<double> wts, int t);
	Torus2D(int N, double L, std::vector<double> sizes, std::vector<double> wts);
	//makes radius=1 disks at pre-set locations
	Torus2D(int N, double L, std::vector<std::vector<double> > locs);

	void set_seed(int seed);
	int get_seed(){return seed;}

	//populates randomly
	void populate();
	//fills with pre-set locations
	void set_populate(std::vector<std::vector<double> > locs);

	void find_set_L(double phi);

	void set_N(int N);
	void set_L(double L);
	void set_1x(int i, double x);
	void set_1y(int i, double y);

	int get_N(){ return N; }
	double get_L(){ return L; }

	double get_1x(int i){return sphs.at(i).get_x();}
	double get_1y(int i){return sphs.at(i).get_y();}

	double get_x(int i){return L*mymath::pmod(sphs.at(i).get_x(),1);}
	double get_y(int i){return L*mymath::pmod(sphs.at(i).get_y(),1);}
	double get_r(int i){return sphs.at(i).get_r();}

	std::vector<Sphere2D> get_sphs(){return sphs;}

	//distance of contact btwn 2 spheres
	double R(int i, int j){return get_r(i)+get_r(j);}

	//gets the kth-quadrant vector from i to j
	double rel_x(int i, int j, int k){
		double relx = mymath::pmod(sphs.at(j).get_x()-sphs.at(i).get_x(),1);
		if(k==1||k==4) return L*relx;
		else return L*(relx-1);
	}
	double rel_y(int i, int j, int k){
		double rely = mymath::pmod(sphs.at(j).get_y()-sphs.at(i).get_y(),1);
		if(k==1||k==2) return L*rely;
		else return L*(rely-1);
	}
	//returns interparticle distance in the k-th quadrant
	double ell(int i, int j, int k);

	//packing fraction
	double pack_frac();

};

#endif /* TORUS2D_H_ */
