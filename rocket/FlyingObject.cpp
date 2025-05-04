#include "FlyingObject.h"


void FlyingObject::computeForces(double vx, double vy, double vz, double x, double y, double z, double& ax, double& ay, double& az, const Enviroment& env)
{
    double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    double rho = env.computeAirDensity(x, y, z);
    double v = sqrt(pow(vx, 2) + pow(vy, 2) + pow(vz, 2));
    ax = (-rho * pow(radius, 2) * pi * shape_coeff * abs(v) * vx + M * U0 * cos(pi / 2) - g * mass * x / r) / mass;
    ay = (-g * mass * y / r - rho * pow(radius, 2) * pi * shape_coeff * abs(v) * vy + M * U0 * sin(pi / 2)) / mass;
    az = (-g * mass * z / r - rho * pow(radius, 2) * pi * shape_coeff * abs(v) * vz + M * U0 * sin(pi / 2)) / mass;


}


FlyingObject::FlyingObject(double mass, double radius, double shape_coeff, double v0, double alpha0, double tetta0, double M, double U0)
    : mass(mass), radius(radius), shape_coeff(shape_coeff), M(M), U0(U0), alpha(alpha0), tetta(tetta0) {
    vx = v0 * cos(alpha0) * sin(tetta0);
    vy = v0 * sin(alpha0) * sin(tetta0);
    vz = v0 * cos(tetta0);
    x = 0;
    y = 0;
    z = 0;
}