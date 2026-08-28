#include <algorithm>
#include <cstdio>
#include <string>

#include "src/aero/LatticeBoltzmannSolver.h"
#include "src/aero/WindTunnel.h"
#include "src/geometry/BoundingBox.h"
#include "src/geometry/Mesh.h"
#include "src/geometry/Transform.h"
#include "src/io/StlLoader.h"
#include "src/utils/Logger.h"

int main(int argc, char** argv) {
    Utils::Logger::Info("Wind Tunnel Modelling Sim -- phase 1 physics harness");

    if (argc < 2) {
        Utils::Logger::Error("Usage: wind_tunnel_physics <path-to-stl> [velocity_mps] [angle_of_attack_deg]");
        return 1;
    }
    std::string stlPath = argv[1];

    Aero::WindTunnel tunnel; // sea-level air, 20 m/s, 2m x 2m x 5m test section
    if (argc >= 3) tunnel.freeStreamVelocity = std::stof(argv[2]);
    float angleOfAttackDeg = (argc >= 4) ? std::stof(argv[3]) : 0.0f;

    Geometry::Mesh mesh;
    std::string err;
    if (!Io::StlLoader::Load(stlPath, mesh, err)) {
        Utils::Logger::Error("Failed to load '" + stlPath + "': " + err);
        return 1;
    }
    Utils::Logger::Info("Loaded mesh: " + std::to_string(mesh.GetTriangleCount()) + " triangles");

    Geometry::BoundingBox localBox = mesh.ComputeBoundingBox();
    if (!localBox.IsValid()) {
        Utils::Logger::Error("Loaded mesh has no geometry");
        return 1;
    }
    Geometry::Vec3 localExtents = localBox.Extents();
    float localMaxExtent = std::max({localExtents.x, localExtents.y, localExtents.z});
    if (localMaxExtent < 1e-9f) {
        Utils::Logger::Error("Loaded mesh is degenerate (zero size)");
        return 1;
    }

    // STL/OBJ carry no units, so auto-fit: scale the model's largest
    // dimension to ~30% of the tunnel's test-section height, then place
    // it a quarter of the way down the tunnel's length (leaving wake
    // room downstream) and centered in Y/Z.
    Geometry::Transform xf;
    xf.scale = (tunnel.testSectionHeight * 0.3f) / localMaxExtent;
    xf.SetAngleOfAttack(angleOfAttackDeg);

    Geometry::Vec3 localCenter = localBox.Center();
    Geometry::Vec3 targetWorldCenter(tunnel.testSectionLength * 0.25f, 0.0f, 0.0f);
    xf.position = targetWorldCenter - xf.TransformDirection(localCenter * xf.scale);

    Aero::LatticeBoltzmannSolver solver;
    solver.Initialize(mesh, xf, tunnel);

    // Note: this is the physically-requested Reynolds number, before the
    // solver's own tau-stability clamp (see the WARN above, if any). The
    // actual simulated Reynolds number is in the final results below.
    Utils::Logger::Info("Stepping solver (requested Re=" +
                         std::to_string(tunnel.ReynoldsNumber(localMaxExtent * xf.scale)) + ")...");

    constexpr int kMaxSteps = 40000;
    constexpr int kProgressInterval = 200;
    int step = 0;
    for (; step < kMaxSteps; ++step) {
        solver.Step();
        if (step % kProgressInterval == 0) {
            Aero::AeroResults r = solver.GetResults();
            Utils::Logger::Info("step " + std::to_string(step) +
                                 ": Cd=" + std::to_string(r.dragCoefficient) +
                                 " Cl=" + std::to_string(r.liftCoefficient) +
                                 " Cm=" + std::to_string(r.momentCoefficient));
        }
        if (solver.HasReachedSteadyState()) {
            Utils::Logger::Info("Converged after " + std::to_string(step) + " steps");
            break;
        }
    }
    bool converged = step < kMaxSteps;
    if (!converged) {
        Utils::Logger::Warn("Hit the " + std::to_string(kMaxSteps) +
                             "-step cap without converging; treat results as provisional.");
    }

    Aero::AeroResults final = solver.GetResults();
    std::printf("\n=== Final Results ===\n");
    std::printf("Steps run:         %d%s\n", step, converged ? "" : " (did not converge)");
    std::printf("Reynolds number:   %.1f\n", final.reynoldsNumber);
    std::printf("Dynamic pressure:  %.3f Pa\n", final.dynamicPressure);
    std::printf("Drag force:        %.4f N   (Cd = %.4f)\n", final.dragForce, final.dragCoefficient);
    std::printf("Lift force:        %.4f N   (Cl = %.4f)\n", final.liftForce, final.liftCoefficient);
    std::printf("Side force:        %.4f N\n", final.sideForce);
    std::printf("Moment coeff:      %.4f (Cm)\n", final.momentCoefficient);
    std::printf("L/D ratio:         %.4f\n", final.LiftToDragRatio());

    return converged ? 0 : 2;
}
