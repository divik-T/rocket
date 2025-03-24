#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "Enviroment.h"
#include "FlyingObject.h"


const double dt = 0.01;
const double OrderOfRK = 4;
const int numOfParam = 6;

/*
class Environment {
public:
    double G, m_e, R_e, rho0, T, R, m_air;

    Environment(double G, double m_e, double R_e, double rho0, double T, double R, double m_air)
        : G(G), m_e(m_e), R_e(R_e), rho0(rho0), T(T), R(R), m_air(m_air) {}

    double computeAirDensity(double y) const {

        return rho0 * exp((-m_air * g * y) / (R * T));
    }
};

*/

/*
class FlyingObject {
public:
    double mass, radius, shape_coeff, x, y, vx, vy, M, U0, alpha;

    FlyingObject(double mass, double radius, double shape_coeff, double v0, double alpha0, double M, double U0)
        : mass(mass), radius(radius), shape_coeff(shape_coeff), M(M), U0(U0), alpha(alpha0) {
        vx = v0 * cos(alpha0);
        vy = v0 * sin(alpha0);
        x = 0;
        y = 0;
    }

    void computeForces(double vx, double vy, double y, double& ax, double& ay, const Environment& env) {
        double rho = env.computeAirDensity(y);


        ax = (-rho * pow(radius, 2) * pi * shape_coeff * vx + M * U0 * cos(alpha)) / mass;
        ay = (-g * mass - rho * pow(radius, 2) * pi * shape_coeff * vy + M * U0 * sin(alpha)) / mass;
    }
};

*/

double ChooseDeltaT(int i, double dt) {
    switch (i) {
    case 0:
        return 0.0;
        break;
    case 1:
        return dt / 2;
        break;
    case 2:
        return dt / 2;
        break;
    case 3:
        return dt;
        break;
    }
}


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
        return 0.0;
        break;
    case 1:
        return 1.0 / 2;
        break;
    case 2:
        return 1.0 / 2;
        break;
    case 3:
        return 1.0;
        break;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void oneStepRungeKutt(
    const std::vector<double>& vecOfLeftParts,
    std::vector<double>& vecOfNewLeftParts,
    const std::vector<double(*)(const std::vector<double>& phaseSpace, const FlyingObject& obj, const Enviroment& env)>& vecOfRightParts,
    const FlyingObject& obj,
    const Enviroment& env,
    double t, double dt) {
    // 
    std::vector <double> stateVector(numOfParam);
    std::vector <double> initialStateVector;
    std::copy(vecOfLeftParts.begin(), vecOfLeftParts.end(), stateVector.begin());
    std::copy(stateVector.begin(), stateVector.end(), vecOfNewLeftParts.begin());
    stateVector[stateVector.size() - 1] = t;
    initialStateVector = stateVector;
    std::vector<double> vecOfDeltaParam(stateVector.size());

    for (int j = 0; j < OrderOfRK; ++j) {
        for (int i = 0; i < vecOfDeltaParam.size() - 1; ++i) {
            vecOfDeltaParam[i] = dt * vecOfRightParts[i](stateVector, obj, env);
        }

        vecOfDeltaParam[vecOfDeltaParam.size() - 1] = dt;

        for (int i = 0; i < vecOfNewLeftParts.size(); ++i) {
            vecOfNewLeftParts[i] += vecOfDeltaParam[i] * bRK(j);
            stateVector[i] = initialStateVector[i] + cRK(j) * vecOfDeltaParam[i];
        }
    }
}

double SecondNewtonLawX(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env) {
    double
        x = stateVector[0],
        y = stateVector[1],
        Vx = stateVector[2],
        Vy = stateVector[3],
        m = stateVector[4],
        t = stateVector[5];

    return (-env.rho * pow(obj.radius, 2) * pi * obj.shape_coeff * Vx + obj.M * obj.U0 * cos(atan(Vy / Vx))) / obj.mass;
}

double SecondNewtonLawY(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env) {
    double
        x = stateVector[0],
        y = stateVector[1],
        Vx = stateVector[2],
        Vy = stateVector[3],
        m = stateVector[4],
        t = stateVector[5];
    return (-g * obj.mass - env.rho * pow(obj.radius, 2) * pi * obj.shape_coeff * Vy + obj.M * obj.U0 * cos(atan(Vy / Vx))) / obj.mass;
}

double VelocityX(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env) {
    double
        x = stateVector[0],
        y = stateVector[1],
        Vx = stateVector[2],
        Vy = stateVector[3],
        m = stateVector[4],
        t = stateVector[5];

    return Vx;
}

double VelocityY(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env) {
    double
        x = stateVector[0],
        y = stateVector[1],
        Vx = stateVector[2],
        Vy = stateVector[3],
        m = stateVector[4],
        t = stateVector[5];

    return Vy;
}

double massPerTime(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env) {
    double
        x = stateVector[0],
        y = stateVector[1],
        Vx = stateVector[2],
        Vy = stateVector[3],
        m = stateVector[4],
        t = stateVector[5];

    double M = obj.M;
    return M;
}

void updateState(std::vector<double>& vecOfNewLeftParts,
    FlyingObject& obj,
    Enviroment& env) {
    obj.x = vecOfNewLeftParts[0];
    obj.y = vecOfNewLeftParts[1];
    obj.vx = vecOfNewLeftParts[2];
    obj.vy = vecOfNewLeftParts[3];
    obj.mass = vecOfNewLeftParts[4];

    obj.alpha = obj.vy / obj.vx;
    env.rho = env.computeAirDensity(obj.y);

}

////////////////////////////////////////////////////////////////////////////////////////////////////

void WriteInFile(FlyingObject& obj, double t, std::ofstream& file) {
    file << t << "\t" << obj.x << "\t" << obj.y << "\n";
    //std::cout << "t: " << t << " | x: " << obj.x << " | y: " << obj.y << std::endl;

}

double TheoreticalTrajectory(double& vy0, double& mass0, double& U0, double& M, double& t) {
    return  vy0 + U0 * log10((mass0) / (mass0 - M * t)) / log10(exp(1)) - g * t;
}

bool СompOfThAndRK(double& vy0, double& mass0, double& U0, double& M, FlyingObject& obj, const Enviroment& env, double t) {
    double v_theor = TheoreticalTrajectory(vy0, mass0, U0, M, t);
    return (abs((obj.vy - v_theor) / (v_theor))) <= 0.01;
}

void CalculationOfTrajectory(FlyingObject& obj, Enviroment& env, double& t, double dt, std::ofstream& file) {
    double z = 0;
    double vy0 = obj.vy;
    double mass0 = obj.mass;
    double Uy0 = obj.U0;
    double M = obj.M;
    bool Correct = 1;
    std::vector<double>newState(numOfParam);
    std::vector<double>State = { obj.x, obj.y, obj.vx, obj.vy, obj.mass };
    std::vector<double(*)(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env)> vecOfFunctions =
    { VelocityX, VelocityY, SecondNewtonLawX, SecondNewtonLawY, massPerTime };

    for (;;) {
        oneStepRungeKutt(State, newState, vecOfFunctions, obj, env, t, dt);
        updateState(newState, obj, env);
        if (obj.y <= 0 || obj.mass <= 0) {
            obj.y = 0;
            obj.x = State[0] + ((State[2] / State[3]) * (0 - State[1]));
            break;
        }
        State = newState;
        if (Correct == 1) {
            Correct = СompOfThAndRK(vy0, mass0, Uy0, M, obj, env, t + dt);
        }
        WriteInFile(obj, t, file);
        if (Correct == 0) {
            СompOfThAndRK(vy0, mass0, Uy0, M, obj, env, t + dt);
        }
        t += dt;
    }
    if (Correct != 1) {
        std::cout << "error > 1% \n";
    }
    std::cout << "Landed in x: " << obj.x << " m in t: " << t << " s\n";
}

int main() {
    const double dt = 0.01;
    Enviroment earth(5.97e24, 6.37e6, 1.23, 237, 8.31, 0.02897);
    FlyingObject stone(5, 0.2, 0, 100, pi / 4, 0.0, 100);

    std::ofstream file("trajectory.txt");
    if (!file) {
        std::cout << "File is not open" << std::endl;
        return 1;
    }

    double t = 0;

    std::cout << "Start  on Earth\n";
    CalculationOfTrajectory(stone, earth, t, dt, file);

    file.close();
    return 0;
}

