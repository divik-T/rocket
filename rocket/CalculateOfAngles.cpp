#include "CalculateOfAngles.h"
#include "Enviroment.h"
#include "FlyingObject.h"
#include "RangeKutt.h"
#include <iostream>

//----------------------------------------------------------------------------------------------------

// ДУБОВСКИЙ

//----------------------------------------------------------------------------------------------------

int numOfParam = 8;
//Начальные собственные уговые скорости в радианах
double p0 = (PI / 180) * 10,
q0 = (PI / 180) * 100,
r0 = (PI / 180) * 10;
// Поворот вектора из локальной в лабораторную систему
Vec3 rotateVectorByQuaternion(const Vec3& v, const Quaternion& q) {
    // Представление вектора как кватерниона с нулевой компонентой w
    Quaternion vectorQuat = { 0.0, v.getX(), v.getY(), v.getZ()};

    // Переворот (конъюгация кватерниона)
    Quaternion qConjugate = { q.getW(), -q.getX(), -q.getY(), -q.getZ()};

    // Применение вращения через кватернион
    Quaternion resultQuat = (q.multiply(vectorQuat)).multiply( qConjugate);

    // Возвращаем результат как новый вектор
    return { resultQuat.getX(), resultQuat.getY(), resultQuat.getZ()};
}


double EulerEqP(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env) {
    double
        p = stateVector[0],
        q = stateVector[1],
        r = stateVector[2];
    double MomentsOfInertia[3];
    
    return (env.getM().getX() / obj.getMomentsOfInertia().getX() - ((obj.getMomentsOfInertia().getZ() - obj.getMomentsOfInertia().getY()) / obj.getMomentsOfInertia().getX()) * q * r);
}

// Правая часть динамического уравнения Эйлера на ось Q
double EulerEqQ(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env) {
    double
        p = stateVector[0],
        q = stateVector[1],
        r = stateVector[2];

    return (env.getM().getY() / obj.getMomentsOfInertia().getY() - ((obj.getMomentsOfInertia().getX() - obj.getMomentsOfInertia().getZ()) / obj.getMomentsOfInertia().getY()) * p * r);
}

// Правая часть динамического уравнения Эйлера на ось R
double EulerEqR(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env) {
    double
        p = stateVector[0],
        q = stateVector[1],
        r = stateVector[2];

    return (env.getM().getZ() / obj.getMomentsOfInertia().getZ() - ((obj.getMomentsOfInertia().getY() - obj.getMomentsOfInertia().getX()) / obj.getMomentsOfInertia().getZ()) * p * q);
}

// Правая часть для угла Phi
double convertLocalToEulerPhi(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env) {
    double p = stateVector[0],
        q = stateVector[1],
        r = stateVector[2],
        phi = stateVector[3],
        psi = stateVector[4],
        teta = stateVector[5];
    double result = (p * sin(psi) + cos(psi) * q) / safeSin(teta);
    return result;
}

// Правая часть для угла Psi
double convertLocalToEulerPsi(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env) {
    double p = stateVector[0],
        q = stateVector[1],
        r = stateVector[2],
        phi = stateVector[3],
        psi = stateVector[4],
        teta = stateVector[5];
    double result = r - (sin(psi) * p + cos(psi) * q) * cos(teta) / safeSin(teta);
    return result;
}

// Правая часть для угла Tetta
double convertLocalToEulerTeta(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env) {
    double p = stateVector[0],
        q = stateVector[1],
        r = stateVector[2],
        phi = stateVector[3],
        psi = stateVector[4],
        teta = stateVector[5];

    double result = (cos(psi) * p - sin(psi) * q);
    return result;
}

void updateState(std::vector<double>& vecOfNewLeftParts,
    FlyingObject& obj,
    Enviroment& env,
    double t, double dt) {
    obj.setAngularVelocities(vecOfNewLeftParts[0], vecOfNewLeftParts[1], vecOfNewLeftParts[2]);
    obj.setEulerAngles(vecOfNewLeftParts[3], vecOfNewLeftParts[5], vecOfNewLeftParts[4]);
    
}

void CalculationOfTrajectory(FlyingObject& obj, Enviroment& env, double& t, double dt, std::ofstream& file) {
    std::vector<double>newState(numOfParam);
    std::vector<double>State = { obj.getAngularVelocities().getX(), obj.getAngularVelocities().getY(), obj.getAngularVelocities().getZ(), obj.getEulerAngles().getX(), obj.getEulerAngles().getZ(), obj.getEulerAngles().getY()};
    std::vector<double(*)(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env)> vecOfFunctions =
    { EulerEqP, EulerEqQ, EulerEqR, convertLocalToEulerPhi, convertLocalToEulerPsi, convertLocalToEulerTeta };
    std::vector<double>Solution(3);
    while (t < 10) {
        oneStepRungeKutt(State, newState, vecOfFunctions, obj, env, t, dt);
        updateState(newState, obj, env, t, dt);
        State = newState;
        t += dt;
        if (std::isnan(obj.getEulerAngles().getX()) || std::isnan(obj.getEulerAngles().getY()) || std::isnan(obj.getEulerAngles().getZ())) {
            std::cerr << "NaN detected at t = " << t << std::endl;
            break;
        }

        WriteInFileRotation(obj, t, file, Solution);
    }
    std::cout << "complete";
}

Vec3 calculateMomentOfSoples(std::map<Vec3, double> coordPercOfOpenSoplesMap, double thrustForceAll) {
    Vec3 resultM = { 0,0,0 };
    for (auto it = coordPercOfOpenSoplesMap.begin(); it != coordPercOfOpenSoplesMap.end(); ++it) {
        const Vec3& radiusVec = it->first;
        double PercOfOpenSople = it->second;
        Vec3 thrustForceDir = { 0,0,1 };
        double thrustForceModule = 0.25 * thrustForceAll * PercOfOpenSople;
        Vec3 thrustForce = thrustForceModule * thrustForceDir;
        Vec3 thrustForceMoment = radiusVec.CrossMult(thrustForce);
        resultM += thrustForceMoment;
    }
    return resultM;
}

std::vector<double> RotationdMotion::calcPercOfOpenSoples() const {
    Vec3 tau;
    tau.setX(x - x0);
    tau.setY(y - y0);
    tau.setZ(z - z0);
    tau.normalize();
    obj.convertVectorLabToLocal(tau);
    std::vector<double> PercOfOpenSoples;
    PercOfOpenSoples.push_back(1 - tau.getX());
    PercOfOpenSoples.push_back(1 + tau.getX());
    PercOfOpenSoples.push_back(1 - tau.getY());
    PercOfOpenSoples.push_back(1 + tau.getY());
    return  PercOfOpenSoples;
}

void WriteInFileRotation(FlyingObject& obj, double t, std::ofstream& file, const std::vector<double>& Solution) {
    file << (obj.getEulerAngles().getX() * 180 / PI) << " " << (obj.getEulerAngles().getY() * 180 / PI) << " " << (obj.getEulerAngles().getZ() * 180 / PI) << std::endl;
}

std::vector<double>analyticalSolution(const FlyingObject& obj, double t) {


    double omega = ((obj.getMomentsOfInertia().getX() - obj.getMomentsOfInertia().getZ()) / obj.getMomentsOfInertia().getX()) * r0;

    // Решение
    double p = q0 * sin(omega * t) + p0 * cos(omega * t);
    double q = q0 * cos(omega * t) - p0 * sin(omega * t);
    double r = r0;

    return { p, q, r };
}