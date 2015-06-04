/*
 * Sphero2D.cpp
 *
 *  Created on: May 27, 2015
 *      Author: izzhov
 */

#include "Sphero2D.h"

Sphero2D::Sphero2D(double x, double y, double u, double v, double r, double a) {
	set_x(x); set_y(y); set_u(u); set_v(v); set_r(r); set_a(a);
}

void Sphero2D::set_x(double x) {this->x = x;}

void Sphero2D::set_y(double y) {this->y = y;}

void Sphero2D::set_u(double u) {this->u = u;}

void Sphero2D::set_v(double v) {this->v = v;}

void Sphero2D::set_r(double r) {this->r = r;}

void Sphero2D::set_a(double a) {this->a = a;}

double Sphero2D::area() {
	return (M_PI * r * r + 4 * (a - 1) * r * r);
}

double Sphero2D::Iom(){
	double A = area(); double m1 = 4*(a-1)*r*r/A; double m2 = M_PI*r*r/A;
	double I1 = (m1/12)*4*r*r*((a-1)*(a-1)+1);
	double I2 = m2*(r*r/2+(a-1)*(a-1)*r*r);
	return I1+I2;
}

void Sphero2D::normalize() {
	double mag = sqrt(u*u+v*v);
	u = u/mag; v = v/mag;
}
