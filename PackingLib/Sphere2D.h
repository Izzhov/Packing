/*
 * Sphere2D.h
 *
 *  Created on: May 19, 2015
 *      Author: root
 */

#ifndef SPHERE2D_H_
#define SPHERE2D_H_

#include <cmath>

class Sphere2D {
private:
	double x;
	double y;
	double r;
public:
	Sphere2D(double x, double y, double r);

	void set_x(double x);
	void set_y(double y);
	void set_r(double r);

	double get_x(){return x;}
	double get_y(){return y;}
	double get_r(){return r;}

	double area();
};



#endif /* SPHERE2D_H_ */
