#include "Logger.h"

#include <chrono>
#include <cstdio>
#include <ctime>

namespace Utils {

namespace {

const char* LevelTag(LogLevel level) {
    switch (level) {
        case LogLevel::Info:    return "INFO";
        case LogLevel::Warning: return "WARN";
        case LogLevel::Error:   return "ERROR";
    }
    return "?";
}

std::string Timestamp() {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tmBuf{};
#if defined(_WIN32)
    localtime_s(&tmBuf, &t);
#else
    localtime_r(&t, &tmBuf);
#endif
    char buf[16];
    std::snprintf(buf, sizeof(buf), "%02d:%02d:%02d", tmBuf.tm_hour, tmBuf.tm_min, tmBuf.tm_sec);
    return std::string(buf);
}

} // namespace

void Logger::Log(LogLevel level, const std::string& message) {
    std::FILE* stream = (level == LogLevel::Info) ? stdout : stderr;
    std::fprintf(stream, "[%s] [%s] %s\n", Timestamp().c_str(), LevelTag(level), message.c_str());
    std::fflush(stream);
}

void Logger::Info(const std::string& message) { Log(LogLevel::Info, message); }
void Logger::Warn(const std::string& message) { Log(LogLevel::Warning, message); }
void Logger::Error(const std::string& message) { Log(LogLevel::Error, message); }

} // namespace Utils
