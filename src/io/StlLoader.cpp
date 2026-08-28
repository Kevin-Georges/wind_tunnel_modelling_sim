#include "StlLoader.h"

#include <cstdint>
#include <fstream>

namespace Io {

using Geometry::Mesh;
using Geometry::Vec3;
using Geometry::Vertex;

bool StlLoader::LoadBinary(const std::string& path, Mesh& outMesh, std::string& err) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        err = "could not open file: " + path;
        return false;
    }

    char header[80];
    file.read(header, sizeof(header));
    if (!file) {
        err = "file too small for a binary STL header";
        return false;
    }

    uint32_t triangleCount = 0;
    file.read(reinterpret_cast<char*>(&triangleCount), sizeof(triangleCount));
    if (!file) {
        err = "failed to read triangle count";
        return false;
    }

    outMesh.vertices.clear();
    outMesh.indices.clear();
    outMesh.vertices.reserve(static_cast<size_t>(triangleCount) * 3);
    outMesh.indices.reserve(static_cast<size_t>(triangleCount) * 3);

    for (uint32_t i = 0; i < triangleCount; ++i) {
        float normal[3];
        float verts[3][3];
        file.read(reinterpret_cast<char*>(normal), sizeof(normal));
        for (auto& v : verts) {
            file.read(reinterpret_cast<char*>(v), sizeof(v));
        }
        uint16_t attributeByteCount = 0;
        file.read(reinterpret_cast<char*>(&attributeByteCount), sizeof(attributeByteCount));
        if (!file) {
            err = "unexpected end of file at triangle " + std::to_string(i) + " of " + std::to_string(triangleCount);
            return false;
        }

        Vec3 faceNormal(normal[0], normal[1], normal[2]);
        uint32_t base = static_cast<uint32_t>(outMesh.vertices.size());
        for (auto& v : verts) {
            outMesh.vertices.push_back(Vertex{Vec3(v[0], v[1], v[2]), faceNormal});
        }
        outMesh.indices.push_back(base);
        outMesh.indices.push_back(base + 1);
        outMesh.indices.push_back(base + 2);
    }

    return true;
}

bool StlLoader::LoadAscii(const std::string& path, Mesh& outMesh, std::string& err) {
    std::ifstream file(path);
    if (!file) {
        err = "could not open file: " + path;
        return false;
    }

    outMesh.vertices.clear();
    outMesh.indices.clear();

    std::string token;
    Vec3 pendingVerts[3];
    int pendingCount = 0;

    while (file >> token) {
        if (token != "vertex") continue;

        float x, y, z;
        if (!(file >> x >> y >> z)) {
            err = "malformed 'vertex' line (expected 3 floats)";
            return false;
        }
        pendingVerts[pendingCount++] = Vec3(x, y, z);

        if (pendingCount == 3) {
            uint32_t base = static_cast<uint32_t>(outMesh.vertices.size());
            for (const auto& p : pendingVerts) {
                outMesh.vertices.push_back(Vertex{p, Vec3()});
            }
            outMesh.indices.push_back(base);
            outMesh.indices.push_back(base + 1);
            outMesh.indices.push_back(base + 2);
            pendingCount = 0;
        }
    }

    if (outMesh.indices.empty()) {
        err = "no triangles found in ASCII STL";
        return false;
    }
    return true;
}

bool StlLoader::Load(const std::string& path, Mesh& outMesh, std::string& err) {
    std::ifstream probe(path, std::ios::binary | std::ios::ate);
    if (!probe) {
        err = "could not open file: " + path;
        return false;
    }
    std::streamsize fileSize = probe.tellg();

    // Binary STL has a fixed-size header + a per-triangle record size, so
    // its total file size is fully determined by the triangle count it
    // declares. Some binary exporters incorrectly start their 80-byte
    // header with the text "solid", so we trust this size check over a
    // "solid" text prefix check.
    bool looksBinary = false;
    if (fileSize >= 84) {
        probe.seekg(80);
        uint32_t count = 0;
        probe.read(reinterpret_cast<char*>(&count), sizeof(count));
        std::streamsize expectedBinarySize = 84 + static_cast<std::streamsize>(count) * 50;
        looksBinary = (expectedBinarySize == fileSize);
    }
    probe.close();

    bool ok = looksBinary ? LoadBinary(path, outMesh, err) : LoadAscii(path, outMesh, err);
    if (!ok && !looksBinary) {
        // Detection guessed ASCII but parsing failed; the file may still
        // be a binary STL that just didn't match the size heuristic.
        std::string fallbackErr;
        if (LoadBinary(path, outMesh, fallbackErr)) {
            ok = true;
        }
    }
    if (!ok) return false;

    outMesh.WeldVertices(1e-5f);
    outMesh.ComputeNormals();
    return true;
}

} // namespace Io
