#include "AnaliticGeometry.h"
#include <cmath>
#include <map>
#include <vector>

double M_PI = 3.1415;

//----------------------------------------------------------------------------------------------------
// Point
//----------------------------------------------------------------------------------------------------

Point::Point(double x, double y, double z)
    : x_(x), y_(y), z_(z) {
}

double Point::distance(const Point& other) const {
    double dx = x_ - other.x_;
    double dy = y_ - other.y_;
    double dz = z_ - other.z_;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

//----------------------------------------------------------------------------------------------------
// Vec3
//----------------------------------------------------------------------------------------------------

Vec3::Vec3(double x, double y, double z)
    : x_(x), y_(y), z_(z) {
}

double Vec3::dot(const Vec3& other) const {
    return x_ * other.x_ + y_ * other.y_ + z_ * other.z_;
}

Vec3 Vec3::CrossMult(const Vec3& b) const {
    return Vec3(
        y_ * b.z_ - z_ * b.y_,
        z_ * b.x_ - x_ * b.z_,
        x_ * b.y_ - y_ * b.x_
    );
}

double Vec3::module() const {
    return sqrt(x_ * x_ + y_ * y_ + z_ * z_);
}

double Vec3::angle(const Vec3& b) const {
    return acos(this->dot(b) / (this->module() * b.module()));
}

void Vec3::normalize() {
    double mod = module();
    if (mod < 1e-10) throw std::runtime_error("Zero-length vector");
    x_ /= mod; y_ /= mod; z_ /= mod;
}

Vec3 Vec3::operator*(double scalar) const {
    return Vec3(x_ * scalar, y_ * scalar, z_ * scalar);
}

Vec3& Vec3::operator+=(const Vec3& other) {
    x_ += other.x_;
    y_ += other.y_;
    z_ += other.z_;
    return *this;
}

bool Vec3::operator<(const Vec3& other) const {
    if (x_ != other.x_) return x_ < other.x_;
    if (y_ != other.y_) return y_ < other.y_;
    return z_ < other.z_;
}

Vec3 operator*(double scalar, const Vec3& v) {
    return v * scalar;
}

//----------------------------------------------------------------------------------------------------
// Quaternion
//----------------------------------------------------------------------------------------------------

Quaternion::Quaternion(double w, double x, double y, double z)
    : w_(w), x_(x), y_(y), z_(z) {
}

void Quaternion::normalize() {
    double length = sqrt(w_ * w_ + x_ * x_ + y_ * y_ + z_ * z_);
    if (length > 0.0) {
        w_ /= length;
        x_ /= length;
        y_ /= length;
        z_ /= length;
    }
}

Quaternion Quaternion::multiply(const Quaternion& q2) const {
    return Quaternion(
        w_ * q2.w_ - x_ * q2.x_ - y_ * q2.y_ - z_ * q2.z_,
        w_ * q2.x_ + x_ * q2.w_ + y_ * q2.z_ - z_ * q2.y_,
        w_ * q2.y_ - x_ * q2.z_ + y_ * q2.w_ + z_ * q2.x_,
        w_ * q2.z_ + x_ * q2.y_ - y_ * q2.x_ + z_ * q2.w_
    );
}

Quaternion Quaternion::conjugate(const Quaternion& q) {
    return Quaternion(q.w_, -q.x_, -q.y_, -q.z_);
}

Quaternion Quaternion::createRotationZ(double angle) {
    double halfAngle = angle / 2.0;
    return Quaternion(cos(halfAngle), 0.0, 0.0, sin(halfAngle));
}

Quaternion Quaternion::createRotationX(double angle) {
    double halfAngle = angle / 2.0;
    return Quaternion(cos(halfAngle), sin(halfAngle), 0.0, 0.0);
}

//----------------------------------------------------------------------------------------------------
// Вращения и преобразования
//----------------------------------------------------------------------------------------------------

Vec3 rotateVectorByQuaternion(const Vec3& v, const Quaternion& q) {
    Quaternion p(0, v.getX(), v.getY(), v.getZ());
    Quaternion result = (q.multiply(p)).multiply(Quaternion::conjugate(q) );
    return Vec3(result.getX(), result.getY(), result.getZ());
}

Vec3 rotateEulerLocalToLab(const Vec3& v, double phi, double theta, double psi) {
    Quaternion qZ = Quaternion::createRotationZ(phi);
    Quaternion qX = Quaternion::createRotationX(theta);
    Quaternion qZ2 = Quaternion::createRotationZ(psi);

    Quaternion qTotal = qZ.multiply( qX.multiply( qZ2));
    qTotal.normalize();

    return rotateVectorByQuaternion(v, qTotal);
}

Vec3 rotateEulerLabToLocal(const Vec3& v, double phi, double theta, double psi) {
    Quaternion qZ = Quaternion::createRotationZ(-psi);
    Quaternion qX = Quaternion::createRotationX(-theta);
    Quaternion qZ2 = Quaternion::createRotationZ(-phi);

    Quaternion qTotal = qZ2.multiply( qX.multiply( qZ));
    qTotal.normalize();

    return rotateVectorByQuaternion(v, qTotal);
}

//----------------------------------------------------------------------------------------------------
// Геометрия цилиндра
//----------------------------------------------------------------------------------------------------

std::map<std::pair<double, double>, Vec3> splitCylinderNormals(double R, double H, double dH, double dGama) {
    std::map<std::pair<double, double>, Vec3> result;
    for (double h = 0; h < H; h += dH) {
        for (double gamma = 0; gamma < 2 * M_PI; gamma += dGama) {
            Vec3 normal(cos(gamma), sin(gamma), 0);
            result[{h, gamma}] = normal;
        }
    }
    return result;
}

Vec3 calculateMomentOfResistance(
    const std::map<std::pair<double, double>, Vec3>& normalMap,
    double resistanceCoeff,
    const Vec3& velocity,
    double R,
    double dS)
{
    Vec3 totalMoment;

    for (std::map<std::pair<double, double>, Vec3>::const_iterator it = normalMap.begin(); it != normalMap.end(); ++it) {
        const std::pair<double, double>& key = it->first;
        const Vec3& normal = it->second;
        Vec3 force = normal * (resistanceCoeff * dS * velocity.module());
        Vec3 radiusVector(R * normal.getX(), R * normal.getY(), 0);
        totalMoment += radiusVector.CrossMult(force);
    }


    return totalMoment;

}


std::map<Vec3, double> createCoordPercOfOpenSoplesMap(double R, const std::vector<double>& PercOfOpenSoples) {
    std::map<Vec3, double> result;
    double angleStep = 2 * M_PI / PercOfOpenSoples.size();
    for (size_t i = 0; i < PercOfOpenSoples.size(); ++i) {
        double angle = i * angleStep;
        Vec3 coord(R * cos(angle), R * sin(angle), 0);
        result[coord] = PercOfOpenSoples[i];
    }
    return result;
}

Vec3 calculateMomentOfSoples(const std::map<Vec3, double>& coordPercOfOpenSoplesMap, double thrustForceAll) {
    Vec3 totalMoment;
    for (const auto& entry : coordPercOfOpenSoplesMap) {
        const Vec3& coord = entry.first;
        double perc = entry.second;
        Vec3 force = coord * (thrustForceAll * perc / 100.0);
        totalMoment += coord.CrossMult(force);
    }
    return totalMoment;
}


double safeSin(double angle) {
    double s = sin(angle);
    return fabs(s) < 1e-6 ? (s >= 0 ? 1e-6 : -1e-6) : s;
}

double normalizePhiPsi(double angle) {
    angle = fmod(angle, 2 * M_PI);
    return angle < 0 ? angle + 2 * M_PI : angle;
}
