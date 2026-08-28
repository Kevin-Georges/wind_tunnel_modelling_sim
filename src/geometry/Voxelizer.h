#pragma once

#include <cstdint>
#include <vector>

#include "BoundingBox.h"
#include "Mesh.h"
#include "Transform.h"
#include "Vec3.h"

namespace Geometry {

struct VoxelGrid {
    int dimX = 0, dimY = 0, dimZ = 0;
    float cellSize = 1.0f;
    Vec3 origin; // world position of the (0,0,0) cell's min corner

    std::vector<uint8_t> occupancy; // 1 = solid, 0 = fluid; flat-indexed

    size_t Index(int x, int y, int z) const {
        return static_cast<size_t>(x) + static_cast<size_t>(dimX) *
               (static_cast<size_t>(y) + static_cast<size_t>(dimY) * static_cast<size_t>(z));
    }

    bool InBounds(int x, int y, int z) const {
        return x >= 0 && y >= 0 && z >= 0 && x < dimX && y < dimY && z < dimZ;
    }

    uint8_t At(int x, int y, int z) const { return occupancy[Index(x, y, z)]; }

    Vec3 CellCenter(int x, int y, int z) const {
        return origin + Vec3((x + 0.5f) * cellSize, (y + 0.5f) * cellSize, (z + 0.5f) * cellSize);
    }
};

class Voxelizer {
public:
    // Voxelizes `mesh` (transformed into world space by `xf`) against a
    // grid spanning `domainBounds` (the tunnel test section) at the given
    // cell size. Marks both the mesh's surface shell and its enclosed
    // interior as solid, so the result is usable directly as a no-slip
    // obstacle mask by the LBM solver.
    static VoxelGrid Voxelize(const Mesh& mesh, const Transform& xf,
                               const BoundingBox& domainBounds, float cellSize);
};

} // namespace Geometry
