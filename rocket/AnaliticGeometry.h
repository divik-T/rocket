#pragma once

#include <cmath>
#include <iostream>
#include <map>
#include <vector>
#include <utility>

//----------------------------------------------------------------------------------------------------
// ДУБОВСКИЙ, ДЕМЯШКЕВИЧ
//----------------------------------------------------------------------------------------------------

class Point {
public:
    Point(double x, double y, double z);

    double distance(const Point& other) const;

    double getX() const { return x_; }
    double getY() const { return y_; }
    double getZ() const { return z_; }

    void setX(double x) { x_ = x; }
    void setY(double y) { y_ = y; }
    void setZ(double z) { z_ = z; }

private:
    double x_, y_, z_;
};

// Класс вектораfdg
class Vec3 {
public:
    Vec3(double x = 0, double y = 0, double z = 0);

    // Геттеры
    double getX() const { return x_; }
    double getY() const { return y_; }
    double getZ() const { return z_; }

    // Сеттеры
    void setX(double x) { x_ = x; }
    void setY(double y) { y_ = y; }
    void setZ(double z) { z_ = z; }

    // Скалярное произведение векторов

    double dot(const Vec3& other) const;

    // Векторное произведение
    Vec3 CrossMult(const Vec3& b) const;

    // Модуль вектора
    double module() const;

    // Расчет угла между векторами
    double angle(const Vec3& b) const;

    // Приведение вектора к единичному
    void normalize();

    // Операторы
    Vec3 operator*(double scalar) const;
    Vec3& operator+=(const Vec3& other);
    bool operator<(const Vec3& other) const;

    // Статический оператор для scalar * Vec3
    friend Vec3 operator*(double scalar, const Vec3& v);

private:
    double x_, y_, z_;
};

class Quaternion {
public:
    Quaternion(double w = 1, double x = 0, double y = 0, double z = 0);

    double getW() const { return w_; }
    double getX() const { return x_; }
    double getY() const { return y_; }
    double getZ() const { return z_; }

    void setW(double w) { w_ = w; }
    void setX(double x) { x_ = x; }
    void setY(double y) { y_ = y; }
    void setZ(double z) { z_ = z; }

    void normalize();

    Quaternion multiply(const Quaternion& q2) const;
    // Статические методы для работы с кватернионами
    
    static Quaternion conjugate(const Quaternion& q);
    static Quaternion createRotationZ(double angle);
    static Quaternion createRotationX(double angle);

private:
    double w_, x_, y_, z_;
};


// Вращение вектора с помощью кватерниона
Vec3 rotateVectorByQuaternion(const Vec3& v, const Quaternion& q);

// Преобразования Эйлера
Vec3 rotateEulerLocalToLab(const Vec3& v, double phi, double theta, double psi);
Vec3 rotateEulerLabToLocal(const Vec3& v, double phi, double theta, double psi);

// Геометрия цилиндра и сопла
std::map<std::pair<double, double>, Vec3> splitCylinderNormals(double R, double H, double dH, double dGama);

Vec3 calculateMomentOfResistance(
    const std::map<std::pair<double, double>, Vec3>& normalMap,
    double resistanceCoeff,
    const Vec3& velocity,
    double R,
    double dS);

std::map<Vec3, double> createCoordPercOfOpenSoplesMap(double R, const std::vector<double>& PercOfOpenSoples);

Vec3 calculateMomentOfSoples(const std::map<Vec3, double>& coordPercOfOpenSoplesMap, double thrustForceAll);

// Вспомогательные функции
double safeSin(double angle);
double normalizePhiPsi(double angle);
