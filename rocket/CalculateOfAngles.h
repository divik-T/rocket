#pragma once
#include "AnaliticGeometry.h"
#include <vector>
#include <fstream>
#include <utility>
#include <map>
#include "FlyingObject.h"


//----------------------------------------------------------------------------------------------------

// дсанбяйхи

//----------------------------------------------------------------------------------------------------


class FlyingObject;

class RotationdMotion {
public:
    double x, x0, y, y0, z, z0;
    const FlyingObject& obj;
    std::vector<double>calcPercOfOpenSoples() const;
};

double EulerEqP(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env);

double EulerEqQ(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env);

double EulerEqR(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env);

double convertLocalToEulerPhi(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env);

double convertLocalToEulerPsi(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env);

double convertLocalToEulerTeta(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env);

void updateState(std::vector<double>&vecOfNewLeftParts,
    FlyingObject & obj,
    Enviroment & env,
    double t, double dt);

void CalculationOfTrajectory(FlyingObject& obj, Enviroment& env, double& t, double dt, std::ofstream& file);


std::vector<double>analyticalSolution(const FlyingObject& obj, double t);

void WriteInFileRotation(FlyingObject& obj, double t, std::ofstream& file, const std::vector<double>& Solution);