#pragma once

#include "Vec3.h"

// Simplified for phase 1: no matrix/quaternion machinery yet since nothing
// renders. Just two named, fixed-axis rotations (yaw about the vertical
// axis, angle of attack about the spanwise axis) plus a position offset --
// enough to place a model in the tunnel and sweep AoA. Revisit with a
// proper ToMatrix()/quaternion representation once rendering (phase 2)
// needs arbitrary orientations.

namespace Geometry {

class Transform {
public:
    Vec3 position;

    // Uniform scale applied before rotation. STL/OBJ files carry no unit
    // information, so an arbitrarily-sized import needs this to fit the
    // tunnel's test section regardless of what units it was modeled in.
    float scale = 1.0f;

    // Applies scale, then rotation (yaw then AoA), then translation.
    Vec3 TransformPoint(const Vec3& localPoint) const;

    // Rotation only (yaw then AoA) -- no scale or translation. Correct
    // for direction/normal vectors, where uniform scale doesn't change
    // direction, only magnitude.
    Vec3 TransformDirection(const Vec3& localDirection) const;

    void SetAngleOfAttack(float degrees) { angleOfAttackDeg_ = degrees; }
    float GetAngleOfAttack() const { return angleOfAttackDeg_; }

    void SetYaw(float degrees) { yawDeg_ = degrees; }
    float GetYaw() const { return yawDeg_; }

private:
    float angleOfAttackDeg_ = 0.0f; // rotation about Z (spanwise)
    float yawDeg_ = 0.0f;           // rotation about Y (vertical)
};

} // namespace Geometry
