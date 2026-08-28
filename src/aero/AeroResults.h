#pragma once

namespace Aero {

struct AeroResults {
    float dragForce = 0.0f;   // N
    float liftForce = 0.0f;   // N
    float sideForce = 0.0f;   // N

    float dragCoefficient = 0.0f;   // Cd
    float liftCoefficient = 0.0f;   // Cl
    float momentCoefficient = 0.0f; // Cm

    float reynoldsNumber = 0.0f;
    float dynamicPressure = 0.0f; // Pa, echoed from WindTunnel at solve time

    int stepIndex = 0;

    float LiftToDragRatio() const {
        if (dragCoefficient == 0.0f) return 0.0f;
        return liftCoefficient / dragCoefficient;
    }
};

} // namespace Aero
