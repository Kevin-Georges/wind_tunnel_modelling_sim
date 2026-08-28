#pragma once

namespace Aero {

struct WindTunnel {
    float freeStreamVelocity = 20.0f;   // m/s
    float angleOfAttackDeg = 0.0f;      // deg
    float yawDeg = 0.0f;                // deg

    float airDensity = 1.225f;          // kg/m^3, sea-level air
    float kinematicViscosity = 1.5e-5f; // m^2/s, sea-level air

    float testSectionWidth = 2.0f;      // m
    float testSectionHeight = 2.0f;     // m
    float testSectionLength = 5.0f;     // m

    float ReynoldsNumber(float characteristicLength) const {
        return freeStreamVelocity * characteristicLength / kinematicViscosity;
    }

    float DynamicPressure() const {
        return 0.5f * airDensity * freeStreamVelocity * freeStreamVelocity;
    }
};

// Model frontal area / test-section cross-sectional area. Real wind
// tunnels start distrusting results somewhere around 5-7.5%; logs a
// warning past that threshold rather than blocking the simulation.
float BlockageRatio(const WindTunnel& tunnel, float frontalArea);
void WarnIfBlocked(const WindTunnel& tunnel, float frontalArea);

} // namespace Aero
