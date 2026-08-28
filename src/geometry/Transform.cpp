#include "Transform.h"

namespace Geometry {

Vec3 Transform::TransformDirection(const Vec3& localDirection) const {
    // Yaw applied first, then angle of attack -- fixed order since the two
    // rotations are about different axes and don't commute.
    Vec3 yawed = RotateY(localDirection, yawDeg_);
    return RotateZ(yawed, angleOfAttackDeg_);
}

Vec3 Transform::TransformPoint(const Vec3& localPoint) const {
    return TransformDirection(localPoint * scale) + position;
}

} // namespace Geometry
