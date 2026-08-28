#include "BoundingBox.h"

#include <algorithm>
#include <limits>

namespace Geometry {

BoundingBox::BoundingBox() {
    float inf = std::numeric_limits<float>::infinity();
    min = Vec3(inf, inf, inf);
    max = Vec3(-inf, -inf, -inf);
    valid = false;
}

void BoundingBox::Union(const Vec3& point) {
    min = Vec3(std::min(min.x, point.x), std::min(min.y, point.y), std::min(min.z, point.z));
    max = Vec3(std::max(max.x, point.x), std::max(max.y, point.y), std::max(max.z, point.z));
    valid = true;
}

void BoundingBox::Union(const BoundingBox& other) {
    if (!other.valid) return;
    Union(other.min);
    Union(other.max);
}

Vec3 BoundingBox::Center() const {
    return (min + max) * 0.5f;
}

Vec3 BoundingBox::Extents() const {
    return max - min;
}

float ComputeFrontalArea(const BoundingBox& box, Axis flowAxis) {
    Vec3 e = box.Extents();
    switch (flowAxis) {
        case Axis::X: return e.y * e.z;
        case Axis::Y: return e.x * e.z;
        case Axis::Z: return e.x * e.y;
    }
    return 0.0f;
}

} // namespace Geometry
