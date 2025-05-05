#pragma once
#include "FlyingObject.h"
class Point {
public:
	double x, y, z;
	Point(double x, double y, double z);
	double distance(const Point& other) const;
};

class Vector {
public:
	double i, j, k;
	Vector(double i , double j, double k);
	Vector(double x1, double x2, double y1, double y2, double z1, double z2);
	double length() const;
	bool ScalarMOfVector( Vector r);
};

class ControlSystem {
public:
	virtual void send(FlyingObject &obj)=0;

};

class ForwardMotion : public ControlSystem {
public:
	double PercentageOfOpening;
	double ForceOfReactivity;
	double MassPerTime0;
	double MassPerTime;
	double U;
	double ComputeForceR(Point initialCords, double r, double phi, double tetta, Point CorrectCords);
	double polynomialApproximation(double deviation);
	void send(FlyingObject& obj) override;
};