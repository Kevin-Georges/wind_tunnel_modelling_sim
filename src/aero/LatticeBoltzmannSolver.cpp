#include "LatticeBoltzmannSolver.h"

#include <algorithm>
#include <cmath>
#include <string>

#include "../geometry/BoundingBox.h"
#include "../utils/Logger.h"

namespace Aero {

using Geometry::Axis;
using Geometry::BoundingBox;
using Geometry::Mesh;
using Geometry::Transform;
using Geometry::Vec3;
using Geometry::VoxelGrid;
using Geometry::Voxelizer;

namespace {

// D3Q19 lattice directions, weights, and opposite-direction table. Index 0
// is the rest particle; 1-6 are the axis directions; 7-18 are the face
// diagonals. Pairs (1,2), (3,4), ... are opposites by construction, which
// is what lets Opposite() below be a one-line formula instead of a table.
constexpr int kDirX[19] = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0};
constexpr int kDirY[19] = {0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1};
constexpr int kDirZ[19] = {0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1};
constexpr float kWeight[19] = {
    1.0f / 3.0f,
    1.0f / 18.0f, 1.0f / 18.0f, 1.0f / 18.0f, 1.0f / 18.0f, 1.0f / 18.0f, 1.0f / 18.0f,
    1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f,
    1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f,
};

int Opposite(int i) {
    if (i == 0) return 0;
    return (i % 2 == 1) ? i + 1 : i - 1;
}

// How many recent steps of drag-coefficient history HasReachedSteadyState()
// looks at. Needs to be long relative to how slowly this solver's forced
// (numerically-stable) viscosity lets the flow field settle -- see the
// comment on tau clamping in Initialize().
constexpr size_t kConvergenceWindow = 800;

} // namespace

float LatticeBoltzmannSolver::Equilibrium(int i, float rho, const Vec3& u) {
    // int * float promotes to float automatically, so this needs no Vec3
    // construction -- this function is the hottest in the solver (~19
    // calls per fluid cell, several times per step), so avoiding a
    // temporary object here measurably matters.
    float eu = kDirX[i] * u.x + kDirY[i] * u.y + kDirZ[i] * u.z;
    float uu = u.Dot(u);
    return kWeight[i] * rho * (1.0f + 3.0f * eu + 4.5f * eu * eu - 1.5f * uu);
}

void LatticeBoltzmannSolver::Initialize(const Mesh& mesh, const Transform& xf, const WindTunnel& tunnel) {
    tunnel_ = tunnel;

    // Flow along +X, test section centered in Y (vertical) and Z (span).
    BoundingBox domain;
    domain.Union(Vec3(0.0f, -tunnel.testSectionHeight * 0.5f, -tunnel.testSectionWidth * 0.5f));
    domain.Union(Vec3(tunnel.testSectionLength, tunnel.testSectionHeight * 0.5f, tunnel.testSectionWidth * 0.5f));

    BoundingBox worldBox;
    for (const auto& v : mesh.vertices) worldBox.Union(xf.TransformPoint(v.position));
    frontalArea_ = std::max(Geometry::ComputeFrontalArea(worldBox, Axis::X), 1e-6f);
    referenceLength_ = std::max(worldBox.Extents().Length(), 1e-6f);
    referencePoint_ = worldBox.Center();

    Aero::WarnIfBlocked(tunnel, frontalArea_);

    // Cell size is derived from the *object's smallest world-space
    // extent*, not a fixed count across the tunnel length. Sizing off
    // tunnel length alone breaks down for elongated/thin objects: a
    // 3:1:1 box auto-fit to the tunnel's height would have its 1-unit
    // cross-section resolved by only ~2-3 cells, producing meaningless
    // drag (confirmed during testing: it read as *higher* drag than a
    // cube of the same longest dimension, the opposite of the real
    // aerodynamic effect of streamlining).
    Vec3 objExtents = worldBox.Extents();
    float smallestFeature = std::max(std::min({objExtents.x, objExtents.y, objExtents.z}), 1e-6f);
    float cellSize = smallestFeature / static_cast<float>(kTargetCellsAcrossSmallestFeature);

    // Cap total cell count for memory/runtime; coarsen (uniformly, so the
    // object's proportions stay correct) rather than silently running out
    // of memory or taking arbitrarily long on a tiny object in a big tunnel.
    Vec3 domainExtents = domain.Extents();
    auto CellCountFor = [&](float cs) -> size_t {
        size_t nx = static_cast<size_t>(std::ceil(domainExtents.x / cs));
        size_t ny = static_cast<size_t>(std::ceil(domainExtents.y / cs));
        size_t nz = static_cast<size_t>(std::ceil(domainExtents.z / cs));
        return nx * ny * nz;
    };
    size_t estimatedCells = CellCountFor(cellSize);
    if (estimatedCells > kMaxCellCount) {
        // Cell count scales as 1/cellSize^3, so scale cellSize by the
        // cube root of the overage to land back under the cap.
        float ratio = static_cast<float>(estimatedCells) / static_cast<float>(kMaxCellCount);
        cellSize *= std::cbrt(ratio);
        Utils::Logger::Warn(
            "Object is small relative to the tunnel; the target resolution would need " +
            std::to_string(estimatedCells) + " cells. Coarsening to keep memory/runtime "
            "bounded (now ~" + std::to_string(CellCountFor(cellSize)) +
            " cells) -- results may be under-resolved.");
    }

    grid_ = Voxelizer::Voxelize(mesh, xf, domain, cellSize);

    // --- Unit conversion: physical (SI) <-> lattice ---
    dx_ = cellSize;
    inletLatticeVelocity_ = 0.05f; // well under the ~0.3 compressibility limit
    dt_ = inletLatticeVelocity_ * dx_ / std::max(tunnel.freeStreamVelocity, 1e-6f);

    float latticeViscosity = tunnel.kinematicViscosity * dt_ / (dx_ * dx_);
    float rawTau = 3.0f * latticeViscosity + 0.5f;

    // Plain BGK collision goes numerically unstable as tau approaches the
    // 0.5 inviscid limit. Real air's Reynolds number for typical model
    // sizes/speeds implies a tau far too close to that limit for this
    // grid resolution to survive -- so clamp to a safe floor and report
    // the (lower) Reynolds number actually being simulated, rather than
    // silently mislabeling results with the unattainable requested Re.
    constexpr float kMinStableTau = 0.6f;
    tau_ = std::max(rawTau, kMinStableTau);
    if (tau_ > rawTau + 1e-6f) {
        effectiveKinematicViscosity_ = (tau_ - 0.5f) / 3.0f * dx_ * dx_ / dt_;
        Utils::Logger::Warn(
            "Requested Reynolds number exceeds what this grid resolution can simulate "
            "stably with BGK collision; clamping tau to " + std::to_string(tau_) +
            " (effective viscosity " + std::to_string(effectiveKinematicViscosity_) +
            " m^2/s vs requested " + std::to_string(tunnel.kinematicViscosity) +
            " m^2/s). Reported Reynolds number reflects what is actually simulated.");
    } else {
        effectiveKinematicViscosity_ = tunnel.kinematicViscosity;
    }

    size_t cellCount = static_cast<size_t>(grid_.dimX) * grid_.dimY * grid_.dimZ;
    f_.assign(cellCount * kQ, 0.0f);
    fNew_.assign(cellCount * kQ, 0.0f);

    // Start the interior at rest (not at free-stream velocity) and let
    // the flow spin up from the inlet. Starting the whole domain --
    // including right around the obstacle -- at full free-stream velocity
    // creates a large, unphysical "impulsive start" transient that takes
    // far longer to relax away than a gentle spin-up from rest does.
    for (int z = 0; z < grid_.dimZ; ++z) {
        for (int y = 0; y < grid_.dimY; ++y) {
            for (int x = 0; x < grid_.dimX; ++x) {
                size_t cell = grid_.Index(x, y, z);
                if (grid_.occupancy[cell]) continue;
                for (int i = 0; i < kQ; ++i) f_[cell * kQ + i] = Equilibrium(i, 1.0f, Vec3());
            }
        }
    }

    stepIndex_ = 0;
    dragHistory_.clear();
    latestResults_ = AeroResults{};
}

void LatticeBoltzmannSolver::Step() {
    const int dimX = grid_.dimX, dimY = grid_.dimY, dimZ = grid_.dimZ;

    // --- Collision: relax every fluid cell's distributions toward local equilibrium ---
    for (int z = 0; z < dimZ; ++z) {
        for (int y = 0; y < dimY; ++y) {
            for (int x = 0; x < dimX; ++x) {
                size_t cell = grid_.Index(x, y, z);
                if (grid_.occupancy[cell]) continue;

                float* fc = &f_[cell * kQ];
                float rho = 0.0f;
                float mx = 0.0f, my = 0.0f, mz = 0.0f;
                for (int i = 0; i < kQ; ++i) {
                    rho += fc[i];
                    mx += kDirX[i] * fc[i];
                    my += kDirY[i] * fc[i];
                    mz += kDirZ[i] * fc[i];
                }
                rho = std::max(rho, 1e-6f);
                Vec3 u(mx / rho, my / rho, mz / rho);

                for (int i = 0; i < kQ; ++i) {
                    float feq = Equilibrium(i, rho, u);
                    fc[i] -= (fc[i] - feq) / tau_;
                }
            }
        }
    }

    // --- Streaming, with boundary handling folded in ---
    std::fill(fNew_.begin(), fNew_.end(), 0.0f);
    pendingForceLattice_ = Vec3();
    pendingTorqueLattice_ = Vec3();

    for (int z = 0; z < dimZ; ++z) {
        for (int y = 0; y < dimY; ++y) {
            for (int x = 0; x < dimX; ++x) {
                size_t cell = grid_.Index(x, y, z);
                if (grid_.occupancy[cell]) continue;

                const float* fc = &f_[cell * kQ];
                for (int i = 0; i < kQ; ++i) {
                    int nx = x + kDirX[i], ny = y + kDirY[i], nz = z + kDirZ[i];
                    bool xOut = nx < 0 || nx >= dimX;
                    bool yOut = ny < 0 || ny >= dimY;
                    bool zOut = nz < 0 || nz >= dimZ;

                    if (xOut) {
                        // Dropped -- the inlet/outlet passes below set x=0
                        // and x=dimX-1 explicitly, overwriting whatever
                        // partial streaming landed there.
                        continue;
                    }

                    if (yOut || zOut) {
                        // Tunnel side/top/bottom wall: no-slip bounce-back.
                        fNew_[cell * kQ + Opposite(i)] += fc[i];
                        continue;
                    }

                    size_t neighbor = grid_.Index(nx, ny, nz);
                    if (grid_.occupancy[neighbor]) {
                        // Obstacle: no-slip bounce-back, and this link
                        // carries momentum exchange with the model.
                        fNew_[cell * kQ + Opposite(i)] += fc[i];

                        Vec3 e(static_cast<float>(kDirX[i]), static_cast<float>(kDirY[i]), static_cast<float>(kDirZ[i]));
                        Vec3 linkForce = e * (2.0f * fc[i]);
                        pendingForceLattice_ += linkForce;
                        Vec3 leverArm = grid_.CellCenter(x, y, z) - referencePoint_;
                        pendingTorqueLattice_ += leverArm.Cross(linkForce);
                    } else {
                        fNew_[neighbor * kQ + i] += fc[i];
                    }
                }
            }
        }
    }

    // Inlet (x=0): Dirichlet, reset to equilibrium at free-stream velocity.
    Vec3 u0(inletLatticeVelocity_, 0.0f, 0.0f);
    for (int z = 0; z < dimZ; ++z) {
        for (int y = 0; y < dimY; ++y) {
            size_t cell = grid_.Index(0, y, z);
            if (grid_.occupancy[cell]) continue;
            for (int i = 0; i < kQ; ++i) fNew_[cell * kQ + i] = Equilibrium(i, 1.0f, u0);
        }
    }

    // Outlet (x=dimX-1): constant-pressure, not raw zero-gradient. A plain
    // copy-from-interior outlet has no way to relieve pressure: it just
    // mirrors whatever built up upstream instead of acting as a mass
    // sink, so the whole domain slowly pressurizes over time and drag
    // decays toward zero as the flow field stagnates (confirmed by
    // instrumenting DebugAverageDensity() during testing -- density rose
    // from 1.0 to 1.14+ and climbing over 10k steps with the old BC).
    // Pinning density to the same reference (1.0) the inlet uses, while
    // taking velocity from the interior neighbor, gives the domain two
    // fixed-pressure anchors so it can reach an actual balance.
    if (dimX >= 2) {
        for (int z = 0; z < dimZ; ++z) {
            for (int y = 0; y < dimY; ++y) {
                size_t outCell = grid_.Index(dimX - 1, y, z);
                if (grid_.occupancy[outCell]) continue;
                size_t inCell = grid_.Index(dimX - 2, y, z);

                const float* fin = &fNew_[inCell * kQ];
                float rhoIn = 0.0f;
                float mx = 0.0f, my = 0.0f, mz = 0.0f;
                for (int i = 0; i < kQ; ++i) {
                    rhoIn += fin[i];
                    mx += kDirX[i] * fin[i];
                    my += kDirY[i] * fin[i];
                    mz += kDirZ[i] * fin[i];
                }
                rhoIn = std::max(rhoIn, 1e-6f);
                Vec3 uIn(mx / rhoIn, my / rhoIn, mz / rhoIn);

                for (int i = 0; i < kQ; ++i) fNew_[outCell * kQ + i] = Equilibrium(i, 1.0f, uIn);
            }
        }
    }

    std::swap(f_, fNew_);
    ++stepIndex_;
    ComputeResults();
}

void LatticeBoltzmannSolver::ComputeResults() {
    // Force scale for converting lattice momentum-exchange to Newtons:
    // rho_physical * dx^4 / dt^2 (mass scale * length scale / time scale^2).
    float forceScale = tunnel_.airDensity * dx_ * dx_ * dx_ * dx_ / (dt_ * dt_);

    Vec3 force = pendingForceLattice_ * forceScale;
    Vec3 torque = pendingTorqueLattice_ * forceScale;

    float q = tunnel_.DynamicPressure();
    float denom = std::max(q * frontalArea_, 1e-9f);

    AeroResults r;
    r.dragForce = force.x;
    r.liftForce = force.y;
    r.sideForce = force.z;
    r.dragCoefficient = force.x / denom;
    r.liftCoefficient = force.y / denom;
    r.momentCoefficient = torque.z / std::max(denom * referenceLength_, 1e-9f);
    r.reynoldsNumber = tunnel_.freeStreamVelocity * referenceLength_ / effectiveKinematicViscosity_;
    r.dynamicPressure = q;
    r.stepIndex = stepIndex_;
    latestResults_ = r;

    dragHistory_.push_back(r.dragCoefficient);
    while (dragHistory_.size() > kConvergenceWindow) dragHistory_.pop_front();
}

float LatticeBoltzmannSolver::DebugAverageDensity() const {
    double sum = 0.0;
    size_t count = 0;
    for (int z = 0; z < grid_.dimZ; ++z) {
        for (int y = 0; y < grid_.dimY; ++y) {
            for (int x = 0; x < grid_.dimX; ++x) {
                size_t cell = grid_.Index(x, y, z);
                if (grid_.occupancy[cell]) continue;
                const float* fc = &f_[cell * kQ];
                float rho = 0.0f;
                for (int i = 0; i < kQ; ++i) rho += fc[i];
                sum += rho;
                ++count;
            }
        }
    }
    return count > 0 ? static_cast<float>(sum / static_cast<double>(count)) : 0.0f;
}

bool LatticeBoltzmannSolver::HasReachedSteadyState() const {
    constexpr int kWarmupSteps = 2000;
    constexpr float kRelativeTolerance = 0.02f;

    if (stepIndex_ < kWarmupSteps) return false;
    if (dragHistory_.size() < kConvergenceWindow) return false;

    float minV = *std::min_element(dragHistory_.begin(), dragHistory_.end());
    float maxV = *std::max_element(dragHistory_.begin(), dragHistory_.end());
    float scale = std::max(std::fabs((minV + maxV) * 0.5f), 1e-6f);
    bool spreadOk = (maxV - minV) / scale < kRelativeTolerance;

    // Spread alone can look small over a momentarily-flat stretch of an
    // otherwise still-declining curve (this genuinely happened during
    // testing: a false "converged" at step ~500 on a curve that kept
    // dropping for another 5000+ steps). Split the window in half and
    // also require the two halves' means to agree, which catches a
    // persistent drift that a pure min/max spread check misses.
    size_t half = dragHistory_.size() / 2;
    float sumFirst = 0.0f, sumSecond = 0.0f;
    for (size_t k = 0; k < dragHistory_.size(); ++k) {
        if (k < half) sumFirst += dragHistory_[k]; else sumSecond += dragHistory_[k];
    }
    float meanFirst = sumFirst / static_cast<float>(half);
    float meanSecond = sumSecond / static_cast<float>(dragHistory_.size() - half);
    bool driftOk = std::fabs(meanSecond - meanFirst) / scale < kRelativeTolerance;

    return spreadOk && driftOk;
}

} // namespace Aero
