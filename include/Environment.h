#ifndef ENVIROMENT_H
#define ENVIROMENT_H

const double g=9.81;
const double pi=3.1415;
const double G=6.67e-11;

class Environment {
    public:
        double  m_e, R_e, rho0, T, R, m_air;
    
        Environment(double m_e, double R_e, double rho0, double T, double R, double m_air);
    
        double computeAirDensity(double y) const;
    };
#endif