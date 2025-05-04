#pragma once
#include <cmath>


const double g = 9.81;
const double pi = acos(-1);
const double G = 6.67e-11;

const double m_e = 5.97e24;
const double R_e = 6.37e6;
const double R = 8.31;

class Enviroment
{
public:

    double    rho0, T, m_air;

    double rho;

    Enviroment(double rho0, double T, double m_air);

    double computeAirDensity(double x, double y, double z) const;


};