
#include <fstream>
#include <cmath>
#include <vector>
#include "Enviroment.h"
#include "FlyingObject.h"
#include "RangeKutt.h"
#include <iostream>

//----------------------------------------------------------------------------------------------------

// ДЕМЯШКЕВИЧ, АНДРЕЙКОВЕЦ, ДИДЕНЬКО

//----------------------------------------------------------------------------------------------------

bool isInteger(double x) {
    return std::fabs(x - std::round(x)) < 1e-3;
}


const double dt = 0.005;

////////////////////////////////////////////////////////////////////////////////////////////////////

//x [0], y [1], z [2], Vx [3], Vy [4], Vz [5], m [6], t [7];

//----------------------------------------------------------------------------------------------------

// АНДРЕЙКОВЕЦ - ПРАВЫЕ ЧАСТИ РУНГЕ-КУТТА

//----------------------------------------------------------------------------------------------------

double SecondNewtonLawX(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env) {
    double
        x = stateVector[0],
        y = stateVector[1],
        z = stateVector[2],
        Vx = stateVector[3],
        Vy = stateVector[4],
        Vz = stateVector[5],
        m = stateVector[6],
        t = stateVector[7];

    return (-env.rho * pow(obj.radius, 2) * pi * obj.shape_coeff * Vx * sqrt(pow(Vx, 2) + pow(Vy, 2) + pow(Vz, 2)) + obj.M * obj.U0 - G * m * m_e * x / pow((sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))), 3)) / m;// изменён угол, потом исправить
}

double SecondNewtonLawY(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env) {
    double
        x = stateVector[0],
        y = stateVector[1],
        z = stateVector[2],
        Vx = stateVector[3],
        Vy = stateVector[4],
        Vz = stateVector[5],
        m = stateVector[6],
        t = stateVector[7];

    return (-G * m * m_e * y / pow((sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))), 3) - env.rho * pow(obj.radius, 2) * pi * obj.shape_coeff * Vy * sqrt(pow(Vx, 2) + pow(Vy, 2) + pow(Vz, 2)) + obj.M * obj.U0 * Vy / sqrt(pow(Vx, 2) + pow(Vy, 2) + pow(Vz, 2))) / m;
}

double SecondNewtonLawZ(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env) {
    double
        x = stateVector[0],
        y = stateVector[1],
        z = stateVector[2],
        Vx = stateVector[3],
        Vy = stateVector[4],
        Vz = stateVector[5],
        m = stateVector[6],
        t = stateVector[7];

    return (-env.rho * pow(obj.radius, 2) * pi * obj.shape_coeff * sqrt(pow(Vx, 2) + pow(Vy, 2) + pow(Vz, 2)) * Vz + obj.M * obj.U0 - G * m * m_e * z / pow((sqrt(pow(x, 2) + pow(y, 2)) + pow(z, 2)), 3)) / m;// изменён угол, потом исправить
}

double VelocityX(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env) {
    double
        x = stateVector[0],
        y = stateVector[1],
        z = stateVector[2],
        Vx = stateVector[3],
        Vy = stateVector[4],
        Vz = stateVector[5],
        m = stateVector[6],
        t = stateVector[7];

    return Vx;
}
//x [0], y [1], z [2], Vx [3], Vy [4], Vz [5], m [6], t [7];
double VelocityY(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env) {
    double
        x = stateVector[0],
        y = stateVector[1],
        z = stateVector[2],
        Vx = stateVector[3],
        Vy = stateVector[4],
        Vz = stateVector[5],
        m = stateVector[6],
        t = stateVector[7];

    return Vy;
}

double VelocityZ(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env) {
    double
        x = stateVector[0],
        y = stateVector[1],
        z = stateVector[2],
        Vx = stateVector[3],
        Vy = stateVector[4],
        Vz = stateVector[5],
        m = stateVector[6],
        t = stateVector[7];

    return Vz;
}

double massPerTime(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env) {
    double
        x = stateVector[0],
        y = stateVector[1],
        z = stateVector[2],
        Vx = stateVector[3],
        Vy = stateVector[4],
        Vz = stateVector[5],
        m = stateVector[6],
        t = stateVector[7];

    double M = obj.M;
    return -M;
}

void updateState(std::vector<double>& vecOfNewLeftParts,
    FlyingObject& obj,
    Enviroment& env) {
    obj.x = vecOfNewLeftParts[0];
    obj.y = vecOfNewLeftParts[1];
    obj.z = vecOfNewLeftParts[2];
    obj.vx = vecOfNewLeftParts[3];
    obj.vy = vecOfNewLeftParts[4];
    obj.vz = vecOfNewLeftParts[5];
    obj.mass = vecOfNewLeftParts[6];

    obj.alpha = obj.vy / obj.vx;
    obj.betta = obj.vz / (sqrt(pow(obj.vx, 2) + pow(obj.vy, 2)));
    env.rho = env.computeAirDensity(obj.x, obj.y, obj.z);

}

////////////////////////////////////////////////////////////////////////////////////////////////////

//----------------------------------------------------------------------------------------------------

// ДЕМЯШКЕВИЧ - ВСЕ НИЖЕ

//----------------------------------------------------------------------------------------------------

void WriteInFile(FlyingObject& obj, double t, std::ofstream& file) {
    file << t << "\t" << obj.x << "\t" << obj.y << "\t" << obj.z << "\n";
    //std::cout << "t: " << t << " | x: " << obj.x << " | y: " << obj.y << std::endl;
}

double TheoreticalTrajectory(double& vy0, double& mass0, double& U0, double& M, double t) {
    return  vy0 + U0 * log10((mass0) / (mass0 - M * t)) / log10(exp(1)) - g * t;
}

bool СompOfThAndRK(double& vy0, double& mass0, double& U0, double& M, FlyingObject& obj, const Enviroment& env, double t) {
    double v_theor = TheoreticalTrajectory(vy0, mass0, U0, M, t);
    return (abs((obj.vy - v_theor) / (v_theor))) <= 0.01;
}

int numOfParam = 8;

void CalculationOfTrajectory(FlyingObject& obj, Enviroment& env, double& t, double dt, std::ofstream& file) {

    double vy0 = obj.vy;
    double mass0 = obj.mass;
    double Uy0 = obj.U0;
    double M = obj.M;
    bool Correct = 1;
    obj.y = 6400000;
    std::vector<double>newState(numOfParam);
    std::vector<double>State = { obj.x, obj.y, obj.z, obj.vx, obj.vy, obj.vz, obj.mass };
    std::vector<double(*)(const std::vector<double>& stateVector, const FlyingObject& obj, const Enviroment& env)> vecOfFunctions =
    { VelocityX, VelocityY, VelocityZ, SecondNewtonLawX, SecondNewtonLawY, SecondNewtonLawZ, massPerTime };

    for (;;) {
        oneStepRungeKutt(State, newState, vecOfFunctions, obj, env, t, dt);
        updateState(newState, obj, env);
        /*if (sqrt(pow(obj.y, 2) + pow(obj.x, 2) + pow(obj.z, 2)) <= 6400000 || obj.mass <= 0) {
            obj.y = 6400000;
            obj.x = State[0] + ((State[2] / State[3]) * (6400000 - State[1]));
            break;
        }*/
        State = newState;

        //std::cout << "runge-kutt Vy: " << obj.vy << " theor Vy: " << TheoreticalTrajectory(vy0, mass0, Uy0, M, t + dt) << " s\n";

        if (Correct == 1) {
            Correct = СompOfThAndRK(vy0, mass0, Uy0, M, obj, env, t + dt);
        }
        if (isInteger(t)) {
            WriteInFile(obj, t, file);
        }
        if (Correct == 0) {
            СompOfThAndRK(vy0, mass0, Uy0, M, obj, env, t + dt);
        }
        t += dt;

    }
    if (Correct != 1) {
        // std::cout << "error > 1% \n";
    }
    std::cout << "Landed in x: " << obj.x << " m in t: " << t << " s\n";
}

int main() {

    const double dt = 0.01;


    const double rho0 = 1.23; // плонтность воздуха у поверзности 
    const double T = 237; // температура
    const double m_air = 0.02897;// масса воздуха 

    Enviroment earth(rho0, T, m_air, 0,0,0);

    const double mass = 5; // масса объекта
    const double radius = 0.2; //радиус камня
    const double shape_coeff = 0; // коэфффициент трения фигуры
    const double v0 = 7000 * sqrt(2); //модуль начальной скорости 
    const double phi = pi / 8; // угол между проекцией начального вектора скорости на плоскость OXY и осью x
    const double betta = pi / 2; //угол между начальным вектором скорости ракеты и осью z (используется система: x направлен вправо, y-вверх, z на нас( перпендикулярно экрану))
    const double M = 0;     //dm/dt
    const double U0 = 1; //скорость истечения газа


    FlyingObject stone(mass, radius, shape_coeff, v0, phi, betta, M, U0, 0,0,0,0,0,0,0,0,0);

    std::ofstream file("../../WindowsProject1/WindowsProject1/trajectory.txt");
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