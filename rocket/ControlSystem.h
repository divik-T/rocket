#pragma once
#include "FlyingObject.h"
#include "AnaliticGeometry.h"

//----------------------------------------------------------------------------------------------------

// ÄÅÌßØÊÅÂÈ×, ÄÓÁÎÂÑÊÈÉ, ÀÍÄÐÅÉÊÎÂÅÖ

//----------------------------------------------------------------------------------------------------
const double MassPerTime0 = 100;

class ControlSystem {
public:
	virtual void send(FlyingObject &obj)=0;

};

class ForwardMotion : public ControlSystem {
	double PercentageOfOpening;
	double ForceOfReactivity;
	double MassPerTime;
	double U;
public:
	double getMassPerTime() const { return MassPerTime; };
	double getForceOfReactivity() const { return ForceOfReactivity; };
	double getU() const { return U; };
	double getPercentageOfOpening() const { return PercentageOfOpening;};

	void setMassPerTime(double  MperT) { this->MassPerTime = MperT; };
	void setForceOfReactivity(double FofR) { this->ForceOfReactivity = FofR; };
	void setU(double  U) { this->U = U; };
	void setPercentageOfOpening(double PofOp) { this->PercentageOfOpening = PofOp; };

	double ComputeForceR(Point initialCords, double r, double phi, double tetta, Point CorrectCords);
	double polynomialApproximation(double deviation);
	void send(FlyingObject& obj) override;
};