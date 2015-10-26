#ifndef __LOGGING_H
#define __LOGGING_H

#define LEVEL_ALL	0
#define LEVEL_TRACE 5
#define LEVEL_DEBUG 10
#define LEVEL_INFO	20
#define LEVEL_ERROR	25
#define LEVEL_OFF	30

#if defined(LOG_LEVEL_ALL)
#define LOG_LEVEL	LEVEL_ALL
#elif defined(LOG_LEVEL_TRACE)
#define LOG_LEVEL	LEVEL_TRACE
#elif defined(LOG_LEVEL_DEBUG)
#define LOG_LEVEL	LEVEL_DEBUG
#elif defined(LOG_LEVEL_INFO)
#define LOG_LEVEL	LEVEL_INFO
#elif defined(LOG_LEVEL_ERROR)
#define LOG_LEVEL	LEVEL_ERROR
#else
#define LOG_LEVEL	LEVEL_OFF
#endif

#if LOG_LEVEL < LEVEL_OFF

#include <iostream>
#define __BASEF__			(strrchr(__FILE__, '\\') ? strrchr(__FILE__, '\\') + 1 : __FILE__)
#define __LOG_PREFIX__		__BASEF__ << ':' << __LINE__ << " - "

#include <iostream>
#define LOG_MESSAGE(level, x)		std::cout << '[' << level << "] " << __LOG_PREFIX__ << x << std::endl;

#if LOG_LEVEL <= LEVEL_TRACE
#define TRACE_ENABLED
#define TRACE_MESSAGE(x)	LOG_MESSAGE("TRACE", x)
#else
#define TRACE_MESSAGE(x)
#endif

#if LOG_LEVEL <= LEVEL_DEBUG
#define DEBUG_ENABLED
#define DEBUG_MESSAGE(x)	LOG_MESSAGE("DEBUG", x)
#else
#define DEBUG_MESSAGE(x)
#endif // LOGLEVEL < DEBUG

#if LOG_LEVEL <= LEVEL_INFO
#define INFO_ENABLED
#define INFO_MESSAGE(x)		LOG_MESSAGE("INFO", x)
#else
#define INFO_MESSAGE(x)
#endif // LOGLEVEL <= INFO

#if LOG_LEVEL <= LEVEL_ERROR
#define ERROR_ENABLED
#define ERROR_MESSAGE(x)	LOG_MESSAGE("ERROR", x)
#else
#define ERROR_MESSAGE(x)
#endif // LOG_LEVEL <= ERROR

#else

#define ERROR_MESSAGE(x)
#define INFO_MESSAGE(x)
#define DEBUG_MESSAGE(x)
#define TRACE_MESSAGE(x)

#endif // LOGLEVEL <= OFF

#endif // __LOGGING_H