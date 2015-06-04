/*
 * Sphero2D.h
 *
 *  Created on: May 27, 2015
 *      Author: izzhov
 */

#ifndef SPHERO2D_H_
#define SPHERO2D_H_

#include <cmath>

class Sphero2D {
private:
	double x;
	double y;
	double u;
	double v;
	double r;
	double a;
public:
	Sphero2D(double x, double y, double u, double v, double r, double a);

	void set_x(double x);
	void set_y(double y);
	void set_u(double u);
	void set_v(double v);
	void set_r(double r);
	void set_a(double a);

	double get_x(){return x;}
	double get_y(){return y;}
	double get_u(){return u;}
	double get_v(){return v;}
	double get_r(){return r;}
	double get_a(){return a;}

	double area();
	double Iom(); //moment of inertia divided by mass, assuming uniform density
	void normalize();
};

#endif /* SPHERO2D_H_ */
