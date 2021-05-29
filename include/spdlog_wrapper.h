#ifndef __SPDLOG_WRAPPER_H__
#define __SPDLOG_WRAPPER_H__
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG  //
#include "spdlog/sinks/rotating_file_sink.h"
#include "spdlog/spdlog.h"

constexpr auto LOGGER_NAME = "log";
constexpr auto LOGGER_SAVE_PATH = "logs/rotating.log";

class FILE_LOGGER {
public:
    FILE_LOGGER()  //: file_logger(spdlog::rotating_logger_mt(LOGGER_NAME, LOGGER_SAVE_PATH, 1048576 * 10, 5))
    {
        spdlog::set_pattern("[%Y/%m/%d %H:%M:%S][%l][%s:%#:%!] %v");
        file_logger = spdlog::rotating_logger_mt(LOGGER_NAME, LOGGER_SAVE_PATH, 1048576 * 10, 5);
        file_logger->set_level(spdlog::level::debug);
    }
    ~FILE_LOGGER(){};
    std::shared_ptr<spdlog::logger> file_logger;
};

class SPDLOGGER {
public:
    static FILE_LOGGER Get_Logger() {
        static FILE_LOGGER log;
        return log;
    }

private:
    SPDLOGGER() = default;
    ~SPDLOGGER() = default;
};

// * Refs https://zhuanlan.zhihu.com/p/337877916
// * Move declearation to above classes, as static object.

#define SPDLOG_DEBUG_FILE(...) SPDLOG_LOGGER_DEBUG(SPDLOGGER::Get_Logger().file_logger, __VA_ARGS__)

#define SPDLOG_INFO_FILE(...) SPDLOG_LOGGER_INFO(SPDLOGGER::Get_Logger().file_logger, __VA_ARGS__)

#define SPDLOG_WARN_FILE(...) SPDLOG_LOGGER_WARN(SPDLOGGER::Get_Logger().file_logger, __VA_ARGS__)

#define SPDLOG_ERROR_FILE(...) SPDLOG_LOGGER_ERROR(SPDLOGGER::Get_Logger().file_logger, __VA_ARGS__)

#define SPDLOG_CRITICAL_FILE(...) SPDLOG_LOGGER_CRITICAL(SPDLOGGER::Get_Logger().file_logger, __VA_ARGS__)

// * Keep for my own reference
// #define SPDLOG_ERROR_FILE(...)                               \
//     auto logger2 = spdlog::get(LOGGER_NAME);                 \
//     if (nullptr == logger2)                                  \
//     {                                                        \
//         auto rotating_logger = spdlog::rotating_logger_mt(   \
//             LOGGER_NAME, LOGGER_SAVE_PATH, 1048576 * 10, 5); \
//         SPDLOG_LOGGER_ERROR(rotating_logger, __VA_ARGS__);   \
//     }                                                        \
//     else                                                     \
//     {                                                        \
//         SPDLOG_LOGGER_ERROR(logger2, __VA_ARGS__);           \
//     }

#endif  //__SPDLOG_WRAPPER_H__
