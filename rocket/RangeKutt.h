#pragma once
#include <vector>
#include "FlyingObject.h"

//----------------------------------------------------------------------------------------------------

// юмдпеийнбеж

//----------------------------------------------------------------------------------------------------

double bRK(int i);

double cRK(int i);

void oneStepRungeKutt(
    const std::vector<double>& vecOfLeftParts,
    std::vector<double>& vecOfNewLeftParts,
    const std::vector<double(*)(const std::vector<double>& phaseSpace, const FlyingObject& obj, const Enviroment& env)>& vecOfRightParts,
    const FlyingObject& obj,
    const Enviroment& env,
    double t, double dt);
