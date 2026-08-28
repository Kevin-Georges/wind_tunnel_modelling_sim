#pragma once

#include <deque>
#include <vector>

#include "../geometry/Mesh.h"
#include "../geometry/Transform.h"
#include "../geometry/Voxelizer.h"
#include "AeroResults.h"
#include "WindTunnel.h"

namespace Aero {

// D3Q19 lattice Boltzmann (BGK) solver. This is the actual physics core:
// it resolves real viscous flow -- separation, wake, pressure and viscous
// drag -- around a voxelized obstacle, unlike an inviscid panel method.
//
// Numerical design, decided up front so the .cpp reads as an
// implementation of a plan rather than a pile of constants:
//
//  - Units: everything internally runs in "lattice units" (dx = one cell,
//    dt = one lattice timestep, rho ~= 1). Initialize() derives dx/dt from
//    the voxel cell size and a target lattice free-stream velocity chosen
//    to stay well under the ~0.3 compressibility limit, then derives the
//    relaxation time tau from the tunnel's physical kinematic viscosity
//    through that dx/dt. Forces are accumulated in lattice units and
//    converted back to Newtons via the standard LBM force scale
//    rho_physical * dx^4 / dt^2 (see Krueger et al., "The Lattice
//    Boltzmann Method", unit conversion chapter).
//
//  - Boundary conditions:
//      * Obstacle (voxelized model) and the tunnel's top/bottom/side
//        walls: no-slip bounce-back.
//      * Inlet (upstream face): Dirichlet -- reset to the equilibrium
//        distribution at the free-stream velocity every step. Simple and
//        stable, if not perfectly non-reflective.
//      * Outlet (downstream face): constant-pressure -- density pinned to
//        the same reference (1.0) the inlet uses, velocity taken from the
//        neighboring interior cell. A plain zero-gradient/copy outlet was
//        tried first and doesn't work: it has no mechanism to relieve
//        pressure, so the domain slowly pressurizes without bound and
//        drag decays toward zero as the flow field stagnates. Pinning
//        both ends to the same reference density is what lets the domain
//        reach an actual mass-balanced steady state.
//
//  - Force/moment extraction: momentum-exchange (bounce-back) method.
//    Each obstacle bounce-back link contributes 2 * e_i * f_i to the
//    force on the body (a population reversing direction transfers twice
//    its momentum), and (leverArm x that force) to the moment about the
//    model's bounding-box center. Only links against the *model*
//    contribute -- tunnel-wall bounce-back is excluded.
//
// Known limitation: single-relaxation-time BGK at a modest grid
// resolution trades accuracy for tractability. Treated as a first
// correctness pass, to be checked against known reference values before
// trusting it on arbitrary geometry.
class LatticeBoltzmannSolver {
public:
    void Initialize(const Geometry::Mesh& mesh, const Geometry::Transform& xf, const WindTunnel& tunnel);

    // Advances exactly one lattice timestep. Takes no dt parameter --
    // LBM's timestep is fixed by the dx/dt chosen in Initialize(), not
    // something a caller can vary per call.
    void Step();

    AeroResults GetResults() const { return latestResults_; }

    bool HasReachedSteadyState() const;

    // Domain-averaged lattice density over fluid cells, for diagnosing
    // mass-conservation drift. Lattice density should stay close to 1.0
    // (what the inlet condition pins it to); a large drift indicates the
    // boundary conditions aren't conserving mass.
    float DebugAverageDensity() const;

private:
    static constexpr int kQ = 19;

    // Cell size is chosen so the object's smallest world-space extent
    // spans this many cells (see Initialize()) -- not a fixed count
    // across the tunnel length, which breaks down for elongated/thin
    // objects. kMaxCellCount caps total grid size for memory/runtime;
    // if the target resolution would exceed it, cell size is coarsened
    // uniformly (with a logged warning) rather than left unbounded.
    static constexpr int kTargetCellsAcrossSmallestFeature = 12;
    static constexpr size_t kMaxCellCount = 1'000'000;

    static float Equilibrium(int i, float rho, const Geometry::Vec3& u);
    void ComputeResults();

    Geometry::VoxelGrid grid_;
    std::vector<float> f_;
    std::vector<float> fNew_;

    WindTunnel tunnel_;
    float dx_ = 1.0f;               // meters per cell
    float dt_ = 1.0f;                // seconds per lattice step
    float tau_ = 1.0f;               // relaxation time, lattice units
    float inletLatticeVelocity_ = 0.05f;

    // The viscosity actually being simulated. Equals the tunnel's real
    // kinematic viscosity unless tau had to be clamped for stability, in
    // which case this is the (higher) effective viscosity that clamp
    // implies -- used to report an honest, achieved Reynolds number
    // instead of the physically-requested but numerically-unstable one.
    float effectiveKinematicViscosity_ = 1.5e-5f;

    float frontalArea_ = 1.0f;       // m^2, for Cd/Cl normalization
    float referenceLength_ = 1.0f;   // m, for Cm normalization and Re
    Geometry::Vec3 referencePoint_;  // world-space moment reference (model bbox center)

    Geometry::Vec3 pendingForceLattice_;
    Geometry::Vec3 pendingTorqueLattice_;

    int stepIndex_ = 0;
    std::deque<float> dragHistory_;
    AeroResults latestResults_;
};

} // namespace Aero
