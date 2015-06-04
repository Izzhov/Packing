/*
 * Sphere2D.cpp
 *
 *  Created on: May 19, 2015
 *      Author: root
 */

#include "Sphere2D.h"

Sphere2D::Sphere2D(double x, double y, double r) {
	set_x(x);
	set_y(y);
	set_r(r);
}

void Sphere2D::set_x(double x) {
	this->x=x;
}

void Sphere2D::set_y(double y) {
	this->y=y;
}

void Sphere2D::set_r(double r) {
	this->r=r;
}

double Sphere2D::area() {
	return M_PI*r*r;
}
