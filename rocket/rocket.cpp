#include <iostream>
#include <fstream>
#include <cmath>
#include "Enviroment.h"
#include "FlyingObject.h"


const double dt = 0.001;


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

void OneStepRungeKutta(FlyingObject& obj, const Enviroment& env, double dt) {

    double k1vx, k2vx, k3vx, k4vx;
    double k1vy, k2vy, k3vy, k4vy;
    double k1x, k2x, k3x, k4x;
    double k1y, k2y, k3y, k4y;
    double ax, ay;

    obj.computeForces(obj.vx, obj.vy, obj.y, ax, ay, env);

    k1vx = dt * ax;
    k1vy = dt * ay;
    k1x = dt * obj.vx;
    k1y = dt * obj.vy;

    obj.computeForces(obj.vx + k1vx / 2, obj.vy + k1vy / 2, obj.y + k1y / 2, ax, ay, env);
    k2vx = dt * ax;
    k2vy = dt * ay;
    k2x = dt * (obj.vx + k1vx / 2);
    k2y = dt * (obj.vy + k1vy / 2);

    obj.computeForces(obj.vx + k2vx / 2, obj.vy + k2vy / 2, obj.y + k2y / 2, ax, ay, env);
    k3vx = dt * ax;
    k3vy = dt * ay;
    k3x = dt * (obj.vx + k2vx / 2);
    k3y = dt * (obj.vy + k2vy / 2);

    obj.computeForces(obj.vx + k3vx, obj.vy + k3vy, obj.y + k3y, ax, ay, env);
    k4vx = dt * ax;
    k4vy = dt * ay;
    k4x = dt * (obj.vx + k3vx);
    k4y = dt * (obj.vy + k3vy);

    obj.vx += (k1vx + 2 * k2vx + 2 * k3vx + k4vx) / 6;
    obj.vy += (k1vy + 2 * k2vy + 2 * k3vy + k4vy) / 6;
    obj.x += (k1x + 2 * k2x + 2 * k3x + k4x) / 6;
    obj.y += (k1y + 2 * k2y + 2 * k3y + k4y) / 6;

    obj.mass = obj.mass - obj.M * dt;
    obj.alpha = atan2(obj.vy, obj.vx);
}

void WriteInFile(FlyingObject& obj, double t, std::ofstream& file) {
    file << t << "\t" << obj.x << "\t" << obj.y << "\n";
    //std::cout << "t: " << t << " | x: " << obj.x << " | y: " << obj.y << std::endl;

}

double TheoreticalTrajectory(double& vy0, double& mass0, double& U0, double& M, double& t) {
    return  vy0 + U0 * log10((mass0) / (mass0 - M * t)) / log10(exp(1)) - g * t;
}

bool СompOfThAndRK(double& vy0, double& mass0, double& U0, double& M, FlyingObject& obj, const Enviroment& env, double& t) {
    double v_theor = TheoreticalTrajectory(vy0, mass0, U0, M, t);
    return (abs((obj.vy - v_theor) / (v_theor))) <= 0.01;
}

void CalculationOfTrajectory(FlyingObject& obj, const Enviroment& env, double& t, double dt, std::ofstream& file) {
    double z = 0;
    double vy0 = obj.vy;
    double mass0 = obj.mass;
    double Uy0 = obj.U0;
    double M = obj.M;
    bool Correct = 1;
    for (;;) {

        OneStepRungeKutta(obj, env, dt);
        if (obj.y <= 0 || obj.mass <= 0) { break; }
        if (Correct == 1) {
            Correct = СompOfThAndRK(vy0, mass0, Uy0, M, obj, env, t);
        }
        WriteInFile(obj, t, file);
        t += dt;
        if (Correct == 0) {
            СompOfThAndRK(vy0, mass0, Uy0, M, obj, env, t);
        }
    }
    if (Correct != 1) {
        std::cout << "error > 1% \n";
    }
    std::cout << "Landed in x: " << obj.x << " m in t: " << t << " s\n";
}



int main() {
    const double dt = 0.01;
    Enviroment earth(5.97e24, 6.37e6, 1.23, 237, 8.31, 0.02897);
    FlyingObject stone(5, 0.2, 0, 100, pi / 4, 0.01, 100);

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
