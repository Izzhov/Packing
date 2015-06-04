/*
 * SpheroPot2D.h
 *
 *  Created on: May 29, 2015
 *      Author: izzhov
 */

#ifndef SPHEROPOT2D_H_
#define SPHEROPOT2D_H_

#include "SpheroTorus2D.h"

#include <vector>
#include <cmath>

class SpheroPot2D {
private:
	double eps; double d;
	double kk; double P;
	SpheroTorus2D& torus;
	double U; std::vector<std::vector<double> > DU; //force
	std::vector<double> Dtau; //torque
	double DL;
	double epstrack;//only for next_eps
public:
	SpheroPot2D(SpheroTorus2D& ptorus);

	double get_U(){return U;}

	//moves and finds new U and F
	void next();

	//calculates energy and force and fills relevant variables
	void calc_U_F();
	void calc_U_F(int i, int j, int k);
	double dUdl(double lijk, double Rij);

	//moves particles
	void move();
	void translate(); //includes volume change
	void rotate();

	//does an iteration while finding eps to check that potential works w/eps
	double next_eps();
	double D_mag_2();//see 6-1-15 paper
	double calc_U();//finds current U *but does not fill*
	double calc_U(int i, int j, int k);
	double get_epstrack(){return epstrack;}
};

#endif /* SPHEROPOT2D_H_ */
