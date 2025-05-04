#include "ControlSystem.h"
#include <cmath>

//Vector
Vector::Vector(double i, double j, double k): i(i), j(j), k(k){}

Vector::Vector(double x1, double x2, double y1, double y2, double z1, double z2): i(x2-x1), j(y2-y1), k(z2-z1){}

double Vector::length() const
{
	return sqrt(i * i + j * j + k * k);
}
//
//Point 
Point::Point(double x, double y, double z):x(x), y(y), z(z){}

double Point::distance(const Point& other) const
{
	return sqrt(pow(x - other.x, 2) + pow(y - other.y, 2) + pow(z - other.z, 2));
}
//
//Control System
bool Vector::ScalarMOfVector( Vector r)
{
	return (this->i*r.i+ this->j*r.j+this->k*r.k)>0;
}

double ControlSystem::polynomialApproximation(double deviation) {
	const double a5 = 0.03978;
	const double a4 = 0.01337;
	const double a3 = -0.02987;
	const double a2 = -0.06837;
	const double a1 = 0.05269;
	const double a0 = 0.4924; 
	return a5 * pow(deviation, 5) + a4 * pow(deviation, 4) +
		a3 * pow(deviation, 3) + a2 * pow(deviation, 2) +
		a1 * deviation + a0;
}


double ControlSystem::ComputeForceR(Point initialCords, double r, double phi, double tetta, Point CorrectCords)
{
	double deviation = initialCords.distance(CorrectCords);
	Vector VToCurrentCords(initialCords.x - CorrectCords.x, initialCords.y - CorrectCords.y, initialCords.z - CorrectCords.z);
	Vector rocket(r * sin(tetta) * cos(phi) - initialCords.x, r * sin(tetta) * sin(phi) - initialCords.y, r * cos(tetta) - initialCords.z);
	deviation /= 1000.0;
	if (rocket.ScalarMOfVector(VToCurrentCords) == 0) {
		deviation = -deviation;
	}
	 MassPerTime =MassPerTime0 * polynomialApproximation(deviation);
	ForceOfReactivity = MassPerTime*U;
	return ForceOfReactivity;

}



