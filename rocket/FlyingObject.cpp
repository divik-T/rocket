#include "FlyingObject.h"


double M_PI = 3.1415;

//----------------------------------------------------------------------------------------------------
// ДИДЕНКО, ДЕМЯШКЕВИЧ, АНДРЕЙКОВЕЦ
//----------------------------------------------------------------------------------------------------

FlyingObject::FlyingObject(double mass, double radius, double shapeCoeff,
    double v0, double alpha0, double betta0,
    double M, double U0,
    double Jp, double Jq, double Jr,
    double phi, double teta, double psi,
    double p, double q, double r)
    : m_mass(mass),
    m_radius(radius),
    m_shapeCoeff(shapeCoeff),
    m_M(M),
    m_U0(U0),
    m_Jp(Jp),
    m_Jq(Jq),
    m_Jr(Jr),
    m_phi(phi),
    m_teta(teta),
    m_psi(psi),
    m_p(p),
    m_q(q),
    m_r(r)
{
    // Инициализация позиции
    setPosition(0, 0, 0);

    // Расчет начальной скорости
    setVelocity(
        v0 * cos(alpha0) * sin(betta0),
        v0 * sin(alpha0) * sin(betta0),
        v0 * cos(betta0)
    );

    // Валидация параметров через сеттеры
    setFlightParams(M, U0, alpha0, betta0);
    setInertiaMoments(Jp, Jq, Jr);
    setEulerAngles(phi, teta, psi);
    setAngularVelocities(p, q, r);
}

//=== Реализация сеттеров ===========================================
void FlyingObject::setPosition(double x, double y, double z) noexcept {
    m_x = x;
    m_y = y;
    m_z = z;
}

void FlyingObject::setVelocity(double vx, double vy, double vz) noexcept {
    m_vx = vx;
    m_vy = vy;
    m_vz = vz;
}

void FlyingObject::setFlightParams(double M, double U0, double alpha, double betta) {
    if (alpha < -M_PI / 2 || alpha > M_PI / 2)
        throw std::invalid_argument("Alpha angle out of range [-π/2, π/2]");
    if (betta < 0 || betta > M_PI)
        throw std::invalid_argument("Betta angle out of range [0, π]");

    m_M = M;
    m_U0 = U0;
    m_alpha = alpha;
    m_betta = betta;
}

void FlyingObject::setInertiaMoments(double Jp, double Jq, double Jr) {
    if (Jp <= 0 || Jq <= 0 || Jr <= 0)
        throw std::invalid_argument("Moments of inertia must be positive");

    m_Jp = Jp;
    m_Jq = Jq;
    m_Jr = Jr;
}

void FlyingObject::setEulerAngles(double phi, double teta, double psi) {
    // Нормализация углов в диапазон [-π, π]
    m_phi = fmod(phi + M_PI, 2 * M_PI) - M_PI;
    m_teta = fmod(teta + M_PI, 2 * M_PI) - M_PI;
    m_psi = fmod(psi + M_PI, 2 * M_PI) - M_PI;
}

void FlyingObject::setAngularVelocities(double p, double q, double r)  {
    m_p = p;
    m_q = q;
    m_r = r;
}

//=== Преобразование систем координат ===============================
void FlyingObject::convertVectorLocalToLab(Vec3& v) const {
    v = rotateEulerLocalToLab(v, m_phi, m_teta, m_psi);
}

void FlyingObject::convertVectorLabToLocal(Vec3& v) const {
    v = rotateEulerLabToLocal(v, m_phi, m_teta, m_psi);
}

//=== Дополнительные методы =========================================
void FlyingObject::applyForces(const Vec3& force) {
    m_vx += force.getX() / m_mass;
    m_vy += force.getY() / m_mass;
    m_vz += force.getZ() / m_mass;
}


void FlyingObject::updateState(double dt) {
    // Интегрирование линейного движения
    m_x += m_vx * dt;
    m_y += m_vy * dt;
    m_z += m_vz * dt;

    // Интегрирование углового движения (упрощенная модель)
    m_phi += m_p * dt;
    m_teta += m_q * dt;
    m_psi += m_r * dt;
}
