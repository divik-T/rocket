#include "ControlSystem.h"
#include <cmath>


double ForwardMotion::polynomialApproximation(double deviation) {
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



double ForwardMotion::ComputeForceR(Point initialCords, double r, double phi, double tetta, Point CorrectCords)
{	
	double deviation = initialCords.distance(CorrectCords);
	Vec3 VToCurrentCords(initialCords.getX() - CorrectCords.getX(), initialCords.getY() - CorrectCords.getY(), initialCords.getZ() - CorrectCords.getZ());
	Vec3 rocket(r * sin(tetta) * cos(phi) - initialCords.getX(), r * sin(tetta) * sin(phi) - initialCords.getY(), r * cos(tetta) - initialCords.getZ());
	deviation /= 1000.0;
	if (rocket.dot(VToCurrentCords) <0) {
		deviation = -deviation;
	}
	MassPerTime =MassPerTime0 * polynomialApproximation(deviation);
	ForceOfReactivity = MassPerTime*U;
	return ForceOfReactivity;
}

void ForwardMotion::send(FlyingObject& obj) {

}
