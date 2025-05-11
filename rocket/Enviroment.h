#pragma once
#include <cmath>
#include "AnaliticGeometry.h"

//----------------------------------------------------------------------------------------------------

// ÄÈÄÅÍÊÎ, ÄÅÌßØÊÅÂÈ×

//----------------------------------------------------------------------------------------------------


const double g = 9.81;
const double pi = acos(-1);
const double PI = pi;

const double G = 6.67e-11;

const double m_e = 5.97e24;
const double R_e = 6.37e6;
const double R = 8.31;
const double rho0 = 1.23; // ïëîíòíîñòü âîçäóõà ó ïîâåðçíîñòè 
const double T0 = 237; // òåìïåðàòóðà
const double m_air = 0.02897;

class Enviroment
{
    double rho, T;
    double  Mp, Mq, Mr;
public:
    double getRho() const { return rho; };
    void setRho(double rho) { this->rho = rho; };

    Vec3 getM() const { return Mp, Mq, Mr; } //????????????!!!!!!!!!!!!!
    void setM(double Mp, double Mq, double Mr) {
        this->Mp = Mp; this->Mr = Mr; this->Mq = Mq;
    }

    double getT() const { return T; };
    void setT(double T) { this->T = T; };

    Enviroment( double Mp, double Mq, double Mr);
    void computeAirDensity(double x, double y, double z) ;

};