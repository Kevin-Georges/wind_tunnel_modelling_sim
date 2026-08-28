#include "Voxelizer.h"

#include <algorithm>
#include <cmath>
#include <queue>

namespace Geometry {

namespace {

int ClampInt(int v, int lo, int hi) { return std::max(lo, std::min(v, hi)); }

// Separating Axis Theorem test: true if `axis` separates the triangle
// (v0,v1,v2, already in box-local space) from the box of `halfSize`
// centered at the origin.
bool AxisSeparates(const Vec3& axis, const Vec3& v0, const Vec3& v1, const Vec3& v2,
                    const Vec3& halfSize) {
    float p0 = v0.Dot(axis);
    float p1 = v1.Dot(axis);
    float p2 = v2.Dot(axis);
    float r = halfSize.x * std::fabs(axis.x) + halfSize.y * std::fabs(axis.y) +
              halfSize.z * std::fabs(axis.z);
    float minP = std::min({p0, p1, p2});
    float maxP = std::max({p0, p1, p2});
    return (minP > r) || (maxP < -r);
}

// Akenine-Moller triangle/AABB overlap test: 3 box-axis tests, 1
// triangle-normal test, 9 edge-cross-box-axis tests. No separating axis
// found among these 13 candidates implies overlap.
bool TriangleIntersectsBox(const Vec3& t0, const Vec3& t1, const Vec3& t2,
                            const Vec3& boxCenter, const Vec3& halfSize) {
    Vec3 v0 = t0 - boxCenter;
    Vec3 v1 = t1 - boxCenter;
    Vec3 v2 = t2 - boxCenter;

    Vec3 edges[3] = {v1 - v0, v2 - v1, v0 - v2};
    Vec3 boxAxes[3] = {Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};

    for (const auto& a : boxAxes) {
        if (AxisSeparates(a, v0, v1, v2, halfSize)) return false;
    }

    for (const auto& e : edges) {
        for (const auto& a : boxAxes) {
            Vec3 axis = e.Cross(a);
            if (axis.LengthSquared() < 1e-12f) continue; // edge parallel to this box axis
            if (AxisSeparates(axis, v0, v1, v2, halfSize)) return false;
        }
    }

    Vec3 normal = edges[0].Cross(edges[1]);
    if (normal.LengthSquared() > 1e-12f) {
        if (AxisSeparates(normal, v0, v1, v2, halfSize)) return false;
    }

    return true;
}

} // namespace

VoxelGrid Voxelizer::Voxelize(const Mesh& mesh, const Transform& xf,
                               const BoundingBox& domainBounds, float cellSize) {
    VoxelGrid grid;
    grid.cellSize = cellSize;
    grid.origin = domainBounds.min;

    Vec3 extents = domainBounds.Extents();
    grid.dimX = std::max(1, static_cast<int>(std::ceil(extents.x / cellSize)));
    grid.dimY = std::max(1, static_cast<int>(std::ceil(extents.y / cellSize)));
    grid.dimZ = std::max(1, static_cast<int>(std::ceil(extents.z / cellSize)));
    grid.occupancy.assign(static_cast<size_t>(grid.dimX) * grid.dimY * grid.dimZ, 0);

    Vec3 halfCell(cellSize * 0.5f, cellSize * 0.5f, cellSize * 0.5f);

    // --- Shell pass: mark cells overlapping the mesh's transformed surface ---
    for (size_t i = 0; i + 2 < mesh.indices.size(); i += 3) {
        Vec3 p0 = xf.TransformPoint(mesh.vertices[mesh.indices[i]].position);
        Vec3 p1 = xf.TransformPoint(mesh.vertices[mesh.indices[i + 1]].position);
        Vec3 p2 = xf.TransformPoint(mesh.vertices[mesh.indices[i + 2]].position);

        BoundingBox triBox;
        triBox.Union(p0);
        triBox.Union(p1);
        triBox.Union(p2);

        int minX = ClampInt(static_cast<int>(std::floor((triBox.min.x - grid.origin.x) / cellSize)), 0, grid.dimX - 1);
        int minY = ClampInt(static_cast<int>(std::floor((triBox.min.y - grid.origin.y) / cellSize)), 0, grid.dimY - 1);
        int minZ = ClampInt(static_cast<int>(std::floor((triBox.min.z - grid.origin.z) / cellSize)), 0, grid.dimZ - 1);
        int maxX = ClampInt(static_cast<int>(std::floor((triBox.max.x - grid.origin.x) / cellSize)), 0, grid.dimX - 1);
        int maxY = ClampInt(static_cast<int>(std::floor((triBox.max.y - grid.origin.y) / cellSize)), 0, grid.dimY - 1);
        int maxZ = ClampInt(static_cast<int>(std::floor((triBox.max.z - grid.origin.z) / cellSize)), 0, grid.dimZ - 1);

        for (int z = minZ; z <= maxZ; ++z) {
            for (int y = minY; y <= maxY; ++y) {
                for (int x = minX; x <= maxX; ++x) {
                    Vec3 cellCenter = grid.CellCenter(x, y, z);
                    if (TriangleIntersectsBox(p0, p1, p2, cellCenter, halfCell)) {
                        grid.occupancy[grid.Index(x, y, z)] = 1;
                    }
                }
            }
        }
    }

    // --- Interior fill: flood-fill from the domain corner (guaranteed
    // exterior to the model, since the model is expected to fit within
    // the domain with margin) to find every cell reachable through fluid
    // without crossing the shell. Anything left unreached and not part
    // of the shell is enclosed interior solid. ---
    std::vector<uint8_t> visited(grid.occupancy.size(), 0);
    std::queue<int> pending;

    auto tryVisit = [&](int x, int y, int z) {
        if (!grid.InBounds(x, y, z)) return;
        size_t idx = grid.Index(x, y, z);
        if (grid.occupancy[idx] == 1) return; // shell cell, don't traverse through it
        if (visited[idx]) return;
        visited[idx] = 1;
        pending.push(static_cast<int>(idx));
    };

    tryVisit(0, 0, 0);
    while (!pending.empty()) {
        int idx = pending.front();
        pending.pop();
        int z = idx / (grid.dimX * grid.dimY);
        int rem = idx % (grid.dimX * grid.dimY);
        int y = rem / grid.dimX;
        int x = rem % grid.dimX;

        tryVisit(x + 1, y, z);
        tryVisit(x - 1, y, z);
        tryVisit(x, y + 1, z);
        tryVisit(x, y - 1, z);
        tryVisit(x, y, z + 1);
        tryVisit(x, y, z - 1);
    }

    for (size_t idx = 0; idx < grid.occupancy.size(); ++idx) {
        if (grid.occupancy[idx] == 1) continue; // already shell
        if (!visited[idx]) grid.occupancy[idx] = 1; // unreached => enclosed interior
    }

    return grid;
}

} // namespace Geometry
