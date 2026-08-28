#pragma once

#include <string>

#include "../geometry/Mesh.h"

namespace Io {

class StlLoader {
public:
    static bool LoadBinary(const std::string& path, Geometry::Mesh& outMesh, std::string& err);
    static bool LoadAscii(const std::string& path, Geometry::Mesh& outMesh, std::string& err);

    // Auto-detects binary vs ASCII, loads, then welds coincident vertices
    // and recomputes smooth vertex normals (STL only carries flat
    // per-facet normals, which can't survive welding as-is).
    static bool Load(const std::string& path, Geometry::Mesh& outMesh, std::string& err);
};

} // namespace Io
