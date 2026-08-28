#pragma once

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "BoundingBox.h"
#include "Vec3.h"

namespace Geometry {

struct Vertex {
    Vec3 position;
    Vec3 normal;
};

class Mesh {
public:
    std::vector<Vertex> vertices;
    std::vector<uint32_t> indices; // triangle list, 3 per face

    // Recomputes per-vertex normals as the normalized sum of adjacent
    // face normals. Overwrites whatever normals are currently stored.
    void ComputeNormals();

    // Merges vertices whose positions are within `epsilon` of each other,
    // remapping indices accordingly. Needed after STL import, which
    // stores an unshared triangle soup.
    void WeldVertices(float epsilon);

    BoundingBox ComputeBoundingBox() const;

    size_t GetTriangleCount() const { return indices.size() / 3; }
    float GetSurfaceArea() const;

    // Divergence-theorem signed volume. Positive for a closed,
    // outward-facing mesh; near zero or negative flags a non-watertight
    // or inverted-winding import worth warning about.
    float GetSignedVolume() const;
};

} // namespace Geometry
