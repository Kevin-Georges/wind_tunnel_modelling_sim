#include "WindTunnel.h"

#include <string>

#include "../utils/Logger.h"

namespace Aero {

float BlockageRatio(const WindTunnel& tunnel, float frontalArea) {
    float crossSection = tunnel.testSectionWidth * tunnel.testSectionHeight;
    if (crossSection <= 0.0f) return 0.0f;
    return frontalArea / crossSection;
}

void WarnIfBlocked(const WindTunnel& tunnel, float frontalArea) {
    constexpr float kBlockageWarningThreshold = 0.075f;
    float ratio = BlockageRatio(tunnel, frontalArea);
    if (ratio > kBlockageWarningThreshold) {
        Utils::Logger::Warn(
            "Blockage ratio " + std::to_string(ratio * 100.0f) +
            "% exceeds the ~7.5% threshold; results may be unreliable. "
            "Consider a larger test section or a smaller model.");
    }
}

} // namespace Aero
