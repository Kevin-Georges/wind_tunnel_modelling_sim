#pragma once

#include "Vec3.h"

namespace Geometry {

enum class Axis { X, Y, Z };

struct BoundingBox {
    Vec3 min;
    Vec3 max;
    bool valid = false;

    BoundingBox();

    void Union(const Vec3& point);
    void Union(const BoundingBox& other);

    Vec3 Center() const;
    Vec3 Extents() const;
    bool IsValid() const { return valid; }
};

// Projected area of the box perpendicular to the tunnel's flow axis.
// A coarse box-silhouette approximation of frontal area -- documented as
// a known simplification versus a true projected-silhouette area.
float ComputeFrontalArea(const BoundingBox& box, Axis flowAxis);

} // namespace Geometry
