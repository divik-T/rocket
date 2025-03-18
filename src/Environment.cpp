#include "../include/Environment.h"
#include<cmath>


Environment::Environment(double m_e, double R_e, double rho0, double T, double R, double m_air)
: m_e(m_e), R_e(R_e), rho0(rho0), T(T), R(R), m_air(m_air) {}

double Environment::computeAirDensity(double y) const {
    return rho0 * exp((-m_air * g * y) / (R * T));
}