#pragma once

#include <string>

namespace Utils {

enum class LogLevel { Info, Warning, Error };

class Logger {
public:
    static void Log(LogLevel level, const std::string& message);
    static void Info(const std::string& message);
    static void Warn(const std::string& message);
    static void Error(const std::string& message);
};

} // namespace Utils
