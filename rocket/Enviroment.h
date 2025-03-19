#pragma once
#include <cmath>


const double g = 9.81;
const double pi = 3.1415;
const double G = 6.67e-11;

class Enviroment
{
public:

    double  m_e, R_e, rho0, T, R, m_air;

    Enviroment(double m_e, double R_e, double rho0, double T, double R, double m_air);

    double computeAirDensity(double y) const;


};

