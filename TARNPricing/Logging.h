#ifndef __LOGGING_H
#define __LOGGING_H

#ifdef LOGGING_ENABLED

#define __BASEF__ (strrchr(__FILE__, '\\') ? strrchr(__FILE__, '\\') + 1 : __FILE__)
#define __LOG_PREFIX__ '[' <<__BASEF__ << ':' << __LINE__ << "] "

#include <iostream>
#define LOG_MESSAGE(x) std::cout << __LOG_PREFIX__ << x << std::endl;
#define LOG_MESSAGE_PARTIAL(x) std::cout << __LOG_PREFIX__ << x;
#define LOG_CR std::cout << std::endl;

#else

#define LOG_MESSAGE(x)
#define LOG_MESSAGE_PARTIAL(x)
#define LOG_CR

#endif // LOGGING_ENABLED

#endif // __LOGGING_H