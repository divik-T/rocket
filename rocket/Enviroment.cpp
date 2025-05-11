#include "Enviroment.h"

//----------------------------------------------------------------------------------------------------

// дхдемйн, делъьйебхв

//----------------------------------------------------------------------------------------------------

void Enviroment::computeAirDensity(double x, double y, double z) 
{
    setRho(  rho0 * exp((-m_air * G * m_e / (pow(R_e, 2)) * sqrt(pow(y, 2) + pow(x, 2) + pow(z, 2))) / (R * T)));
}




Enviroment::Enviroment(double Mp, double Mq, double Mr)
    :  rho(rho0), Mp(Mp), Mq(Mq), Mr(Mr) {
}