#ifndef FLYINGOBJECT_H
#define FLYINGOBJECT_H

#include <cmath>
#include "Environment.h"

class FlyingObject {
    public:
        double mass, radius, shape_coeff, x, y, vx, vy, M, U0, alpha;
    
        FlyingObject(double mass, double radius, double shape_coeff, double v0, double alpha0, double M, double U0);
    
        void computeForces(double vx, double vy, double y, double& ax, double& ay, const Environment& env);
    };
#endif