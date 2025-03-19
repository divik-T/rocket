#include "FlyingObject.h"


void FlyingObject::computeForces(double vx, double vy, double y, double& ax, double& ay, const Enviroment& env)
{
    double rho = env.computeAirDensity(y);
    ax = (-rho * pow(radius, 2) * pi * shape_coeff * abs(vx) * vx + M * U0 * cos(pi / 2)) / mass;
    ay = (-g * mass - rho * pow(radius, 2) * pi * shape_coeff * abs(vy) * vy + M * U0 * sin(pi / 2)) / mass;


}


FlyingObject::FlyingObject(double mass, double radius, double shape_coeff, double v0, double alpha0, double M, double U0)
    : mass(mass), radius(radius), shape_coeff(shape_coeff), M(M), U0(U0), alpha(alpha0) {
    vx = v0 * cos(alpha0);
    vy = v0 * sin(alpha0);
    x = 0;
    y = 0;
}
