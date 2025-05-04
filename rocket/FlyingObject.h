#pragma once
#include "Enviroment.h"


class FlyingObject
{
public:

    double mass, radius, shape_coeff, x, y, z, vx, vy, vz, M, U0, alpha, tetta;

    FlyingObject(double mass, double radius, double shape_coeff, double v0, double alpha0, double tetta0, double M, double U0);

    void computeForces(double vx, double vy, double vz, double x, double y, double z, double& ax, double& ay, double& az, const Enviroment& env);


};
