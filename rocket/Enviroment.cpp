#include "Enviroment.h"

double Enviroment::computeAirDensity(double x, double y, double z) const
{
    return rho0 * exp((-m_air * g * sqrt(pow(y, 2) + pow(x, 2) + pow(z, 2))) / (R * T));
}



Enviroment::Enviroment(double rho0, double T, double m_air)
    : rho0(rho0), T(T), m_air(m_air), rho(rho0) {
}