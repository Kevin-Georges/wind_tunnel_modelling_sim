#pragma once

#include <cmath>

// Minimal vector math for phase 1. GLM was scaffolded originally, but the
// physics core only needs vec3 arithmetic and a couple of fixed-axis
// rotations -- pulling in a dependency for that would be the kind of
// nice-to-have we're deliberately deferring. Revisit if/when phase 2
// (rendering, arbitrary-axis transforms) needs matrices/quaternions.

namespace Geometry {

struct Vec3 {
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;

    Vec3() = default;
    Vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator+(const Vec3& o) const { return Vec3(x + o.x, y + o.y, z + o.z); }
    Vec3 operator-(const Vec3& o) const { return Vec3(x - o.x, y - o.y, z - o.z); }
    Vec3 operator-() const { return Vec3(-x, -y, -z); }
    Vec3 operator*(float s) const { return Vec3(x * s, y * s, z * s); }
    Vec3 operator/(float s) const { return Vec3(x / s, y / s, z / s); }

    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }

    float Dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }

    Vec3 Cross(const Vec3& o) const {
        return Vec3(y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x);
    }

    float LengthSquared() const { return x * x + y * y + z * z; }
    float Length() const { return std::sqrt(LengthSquared()); }

    Vec3 Normalized() const {
        float len = Length();
        if (len < 1e-12f) return Vec3(0.0f, 0.0f, 0.0f);
        return *this / len;
    }
};

inline Vec3 operator*(float s, const Vec3& v) { return v * s; }

constexpr float kPi = 3.14159265358979323846f;
inline float DegToRad(float degrees) { return degrees * (kPi / 180.0f); }

// Rotation about the world Y (vertical) axis -- used for yaw.
inline Vec3 RotateY(const Vec3& v, float degrees) {
    float rad = DegToRad(degrees);
    float c = std::cos(rad), s = std::sin(rad);
    return Vec3(c * v.x + s * v.z, v.y, -s * v.x + c * v.z);
}

// Rotation about the world Z (spanwise) axis -- used for angle of attack.
inline Vec3 RotateZ(const Vec3& v, float degrees) {
    float rad = DegToRad(degrees);
    float c = std::cos(rad), s = std::sin(rad);
    return Vec3(c * v.x - s * v.y, s * v.x + c * v.y, v.z);
}

} // namespace Geometry
