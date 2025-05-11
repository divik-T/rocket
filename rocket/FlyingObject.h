#pragma once
#include "Enviroment.h"
#include "AnaliticGeometry.h"
#include "CalculateOfAngles.h"

//----------------------------------------------------------------------------------------------------
// ДЕМЯШКЕВИЧ, ДУБОВСКИЙ
//----------------------------------------------------------------------------------------------------

class FlyingObject {
private:
    // Основные физические параметры
    double m_mass;
    double m_radius;
    double m_shapeCoeff;

    // Положение в пространстве
    double m_x;
    double m_y;
    double m_z;

    // Линейные скорости
    double m_vx;
    double m_vy;
    double m_vz;

    // Параметры движения
    double m_M;
    double m_U0;
    double m_alpha;
    double m_betta;

    // Моменты инерции
    double m_Jp;
    double m_Jq;
    double m_Jr;

    // Углы Эйлера
    double m_phi;    // Угол прецессии
    double m_teta;   // Угол нутации
    double m_psi;    // Угол собственного вращения

    // Угловые скорости
    double m_p;
    double m_q;
    double m_r;

public:
    FlyingObject(double mass, double radius, double shapeCoeff,
        double v0, double alpha0, double betta0,
        double M, double U0,
        double Jp, double Jq, double Jr,
        double phi, double teta, double psi,
        double p, double q, double r);

    //=== Геттеры ==================================================
    double getMass() const noexcept { return m_mass; }
    double getRadius() const noexcept { return m_radius; }
    double getShapeCoeff() const noexcept { return m_shapeCoeff; }

    Vec3 getPosition() const noexcept { return { m_x, m_y, m_z }; }
    Vec3 getVelocity() const noexcept { return { m_vx, m_vy, m_vz }; }

    double getM() const noexcept { return m_M; }
    double getU0() const noexcept { return m_U0; }
    double getAlpha() const noexcept { return m_alpha; }
    double getBetta() const noexcept { return m_betta; }

    Vec3 getMomentsOfInertia() const noexcept { return { m_Jp, m_Jq, m_Jr }; }
    Vec3 getEulerAngles() const noexcept { return { m_phi, m_teta, m_psi }; }
    Vec3 getAngularVelocities() const noexcept { return { m_p, m_q, m_r }; }

    //=== Сеттеры с базовой валидацией =============================
    void setPosition(double x, double y, double z) noexcept;
    void setVelocity(double vx, double vy, double vz) noexcept;
    void setFlightParams(double M, double U0, double alpha, double betta);
    void setInertiaMoments(double Jp, double Jq, double Jr);
    void setEulerAngles(double phi, double teta, double psi);
    void setAngularVelocities(double p, double q, double r);

    //=== Методы преобразования координат ==========================
    void convertVectorLocalToLab(Vec3& v) const;
    void convertVectorLabToLocal(Vec3& v) const;

    //=== Системные методы =========================================
    void updateState(double dt);  // Для интеграции состояний
    void applyForces(const Vec3& force);
    void applyTorques(const Vec3& torque);
};
