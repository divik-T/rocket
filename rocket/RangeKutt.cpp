#include "RangeKutt.h"
#include <iostream>


//----------------------------------------------------------------------------------------------------

// ¿Õƒ–≈… Œ¬≈÷

//----------------------------------------------------------------------------------------------------


int OrderOfRK = 4;
int numOfParam = 7;

double bRK(int i) {
    switch (i) {
    case 0:
        return 1.0 / 6;
        break;
    case 1:
        return 2.0 / 6;
        break;
    case 2:
        return 2.0 / 6;
        break;
    case 3:
        return 1.0 / 6;
        break;
    }
}

double cRK(int i) {
    switch (i) {
    case 0:
        return 0.5;
        break;
    case 1:
        return 0.5;
        break;
    case 2:
        return 1;
        break;
    case 3:
        return 1.0;
        break;
    }
}

void oneStepRungeKutt(
    const std::vector<double>& vecOfLeftParts,
    std::vector<double>& vecOfNewLeftParts,
    const std::vector<double(*)(const std::vector<double>& phaseSpace, const FlyingObject& obj, const Enviroment& env)>& vecOfRightParts,
    const FlyingObject& obj,
    const Enviroment& env,
    double t, double dt) {

    std::vector <double> stateVector(vecOfNewLeftParts.size());
    std::vector <double> initialStateVector;
    std::copy(vecOfLeftParts.begin(), vecOfLeftParts.end(), stateVector.begin());
    std::copy(stateVector.begin(), stateVector.end(), vecOfNewLeftParts.begin());
    stateVector[stateVector.size() - 1] = t;
    initialStateVector = stateVector;
    std::vector<double> vecOfDeltaParam(stateVector.size());
    vecOfDeltaParam[vecOfDeltaParam.size() - 1] = dt;

    for (int j = 0; j < OrderOfRK; ++j) {
        for (int i = 0; i < vecOfDeltaParam.size() - 1; ++i) {
            vecOfDeltaParam[i] = dt * vecOfRightParts[i](stateVector, obj, env);
        }

        for (int i = 0; i < vecOfNewLeftParts.size(); ++i) {
            vecOfNewLeftParts[i] += vecOfDeltaParam[i] * bRK(j);
            stateVector[i] = initialStateVector[i] + cRK(j) * vecOfDeltaParam[i];
        }
    }
}
