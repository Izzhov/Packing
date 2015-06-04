/*
 * Potential2DNoNbr.h
 *
 *  Created on: May 26, 2015
 *      Author: izzhov
 */

#ifndef POTENTIAL2DNONBR_H_
#define POTENTIAL2DNONBR_H_

#include "Torus2D.h"

#include <vector>
#include <cmath>

class Potential2DNoNbr {
private:
	double eps; double d;
	double kk; double P;
	Torus2D& torus;
	double U; std::vector<std::vector<double> > DU; double DL;
public:
	Potential2DNoNbr(Torus2D& torus);

	double get_U(){return U;}
	std::vector<std::vector<double> > get_DU(){return DU;}
	double get_DL(){return DL;}

	double D_mag_2();

	void harm_osc();

	void grad_harm_osc();
	double dxj_ell(int i, int j, int k); // d/dxj(ellijk)
	double dxi_ell(int i, int j, int k);
	double dyj_ell(int i, int j, int k);
	double dyi_ell(int i, int j, int k);

	void L_harm_osc();
	double dL_ell(int i, int j, int k);

	void full_next();

	//stores new coords and U, but not new DU or DL
	void full_next_debug();
	double next_epscheck();
	void next_rest();
	void move();

	double get_P(){return P;}

	void set_L(double L){torus.set_L(L);} //for debugging object

	bool is_done(); //check if the total force is < d of pressure
	//and actual pressure is < twice & > half assigned press
	//see notebook (Press) for deets
};

#endif /* POTENTIAL2DNONBR_H_ */
