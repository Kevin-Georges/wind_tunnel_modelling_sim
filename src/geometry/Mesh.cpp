#include "Mesh.h"

#include <cmath>

namespace Geometry {

void Mesh::ComputeNormals() {
    for (auto& v : vertices) {
        v.normal = Vec3(0.0f, 0.0f, 0.0f);
    }

    for (size_t i = 0; i + 2 < indices.size(); i += 3) {
        uint32_t i0 = indices[i], i1 = indices[i + 1], i2 = indices[i + 2];
        const Vec3& p0 = vertices[i0].position;
        const Vec3& p1 = vertices[i1].position;
        const Vec3& p2 = vertices[i2].position;

        // Unnormalized face normal; its magnitude is proportional to
        // triangle area, so larger triangles naturally weight the
        // average more -- a common, cheap approximation of area-weighted
        // vertex normals.
        Vec3 faceNormal = (p1 - p0).Cross(p2 - p0);

        vertices[i0].normal += faceNormal;
        vertices[i1].normal += faceNormal;
        vertices[i2].normal += faceNormal;
    }

    for (auto& v : vertices) {
        v.normal = v.normal.Normalized();
    }
}

namespace {

struct CellKey {
    int32_t x, y, z;
    bool operator==(const CellKey& o) const { return x == o.x && y == o.y && z == o.z; }
};

struct CellKeyHash {
    size_t operator()(const CellKey& k) const {
        // Simple mix; collisions just mean an extra distance check, not
        // a correctness problem.
        size_t h = static_cast<size_t>(k.x) * 73856093u;
        h ^= static_cast<size_t>(k.y) * 19349663u;
        h ^= static_cast<size_t>(k.z) * 83492791u;
        return h;
    }
};

CellKey KeyFor(const Vec3& p, float epsilon) {
    return CellKey{
        static_cast<int32_t>(std::floor(p.x / epsilon)),
        static_cast<int32_t>(std::floor(p.y / epsilon)),
        static_cast<int32_t>(std::floor(p.z / epsilon)),
    };
}

} // namespace

void Mesh::WeldVertices(float epsilon) {
    if (epsilon <= 0.0f) return;

    std::unordered_map<CellKey, std::vector<uint32_t>, CellKeyHash> buckets;
    std::vector<Vertex> newVertices;
    newVertices.reserve(vertices.size());
    std::vector<uint32_t> remap(vertices.size());

    float epsilonSq = epsilon * epsilon;

    for (size_t oldIdx = 0; oldIdx < vertices.size(); ++oldIdx) {
        const Vertex& v = vertices[oldIdx];
        CellKey key = KeyFor(v.position, epsilon);

        // Only checks the exact cell a point falls in, not neighboring
        // cells -- a near-duplicate that lands just across a cell
        // boundary won't be merged. Acceptable approximation for phase 1;
        // revisit with a 27-neighbor search if STL imports show visible
        // seams from unmerged boundary vertices.
        uint32_t foundIdx = UINT32_MAX;
        auto it = buckets.find(key);
        if (it != buckets.end()) {
            for (uint32_t candidate : it->second) {
                Vec3 diff = newVertices[candidate].position - v.position;
                if (diff.LengthSquared() <= epsilonSq) {
                    foundIdx = candidate;
                    break;
                }
            }
        }

        if (foundIdx == UINT32_MAX) {
            foundIdx = static_cast<uint32_t>(newVertices.size());
            newVertices.push_back(v);
            buckets[key].push_back(foundIdx);
        }

        remap[oldIdx] = foundIdx;
    }

    for (auto& idx : indices) {
        idx = remap[idx];
    }
    vertices = std::move(newVertices);
}

BoundingBox Mesh::ComputeBoundingBox() const {
    BoundingBox box;
    for (const auto& v : vertices) {
        box.Union(v.position);
    }
    return box;
}

float Mesh::GetSurfaceArea() const {
    float area = 0.0f;
    for (size_t i = 0; i + 2 < indices.size(); i += 3) {
        const Vec3& p0 = vertices[indices[i]].position;
        const Vec3& p1 = vertices[indices[i + 1]].position;
        const Vec3& p2 = vertices[indices[i + 2]].position;
        area += 0.5f * (p1 - p0).Cross(p2 - p0).Length();
    }
    return area;
}

float Mesh::GetSignedVolume() const {
    // Sum of signed tetrahedra formed by each triangle and the origin
    // (divergence theorem). Requires consistently outward-facing winding.
    float volume = 0.0f;
    for (size_t i = 0; i + 2 < indices.size(); i += 3) {
        const Vec3& p0 = vertices[indices[i]].position;
        const Vec3& p1 = vertices[indices[i + 1]].position;
        const Vec3& p2 = vertices[indices[i + 2]].position;
        volume += p0.Dot(p1.Cross(p2)) / 6.0f;
    }
    return volume;
}

} // namespace Geometry
