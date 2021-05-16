#ifndef __SPDLOG_WRAPPER_H__
#define __SPDLOG_WRAPPER_H__
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_INFO //
#include "spdlog/spdlog.h"
#include "spdlog/sinks/rotating_file_sink.h"

constexpr auto LOGGER_NAME = "log";
constexpr auto LOGGER_SAVE_PATH = "logs/rotating.log";

// * Copy from https://zhuanlan.zhihu.com/p/337877916

#define SPDLOG_DEBUG_FILE(...)                                   \
    {                                                            \
        auto logger0 = spdlog::get(LOGGER_NAME);                 \
        if (nullptr == logger0)                                  \
        {                                                        \
            auto rotating_logger = spdlog::rotating_logger_mt(   \
                LOGGER_NAME, LOGGER_SAVE_PATH, 1048576 * 10, 5); \
            SPDLOG_LOGGER_DEBUG(rotating_logger, __VA_ARGS__);   \
        }                                                        \
        else                                                     \
        {                                                        \
            SPDLOG_LOGGER_DEBUG(logger0, __VA_ARGS__);           \
        }                                                        \
    }

#define SPDLOG_INFO_FILE(...)                                    \
    {                                                            \
        auto logger0 = spdlog::get(LOGGER_NAME);                 \
        if (nullptr == logger0)                                  \
        {                                                        \
            auto rotating_logger = spdlog::rotating_logger_mt(   \
                LOGGER_NAME, LOGGER_SAVE_PATH, 1048576 * 10, 5); \
            SPDLOG_LOGGER_INFO(rotating_logger, __VA_ARGS__);    \
        }                                                        \
        else                                                     \
        {                                                        \
            SPDLOG_LOGGER_INFO(logger0, __VA_ARGS__);            \
        }                                                        \
    }

#define SPDLOG_WARN_FILE(...)                                    \
    {                                                            \
        auto logger1 = spdlog::get(LOGGER_NAME);                 \
        if (nullptr == logger1)                                  \
        {                                                        \
            auto rotating_logger = spdlog::rotating_logger_mt(   \
                LOGGER_NAME, LOGGER_SAVE_PATH, 1048576 * 10, 5); \
            SPDLOG_LOGGER_WARN(rotating_logger, __VA_ARGS__);    \
        }                                                        \
        else                                                     \
        {                                                        \
            SPDLOG_LOGGER_WARN(logger1, __VA_ARGS__);            \
        }                                                        \
    }

#define SPDLOG_ERROR_FILE(...)                               \
    auto logger2 = spdlog::get(LOGGER_NAME);                 \
    if (nullptr == logger2)                                  \
    {                                                        \
        auto rotating_logger = spdlog::rotating_logger_mt(   \
            LOGGER_NAME, LOGGER_SAVE_PATH, 1048576 * 10, 5); \
        SPDLOG_LOGGER_ERROR(rotating_logger, __VA_ARGS__);   \
    }                                                        \
    else                                                     \
    {                                                        \
        SPDLOG_LOGGER_ERROR(logger2, __VA_ARGS__);           \
    }

#define SPDLOG_CRITICAL_FILE(...)                             \
    auto logger3 = spdlog::get(LOGGER_NAME);                  \
    if (nullptr == logger3)                                   \
    {                                                         \
        auto rotating_logger = spdlog::rotating_logger_mt(    \
            LOGGER_NAME, LOGGER_SAVE_PATH, 1048576 * 10, 5);  \
        SPDLOG_LOGGER_CRITICAL(rotating_logger, __VA_ARGS__); \
    }                                                         \
    else                                                      \
    {                                                         \
        SPDLOG_LOGGER_CRITICAL(logger3, __VA_ARGS__);         \
    }

#endif //__SPDLOG_WRAPPER_H__
