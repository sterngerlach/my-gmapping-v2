
/* memory_usage.hpp */

#ifndef MY_GMAPPING_MEMORY_USAGE_HPP
#define MY_GMAPPING_MEMORY_USAGE_HPP

#include <cstdint>

namespace MyGMapping {

/* Get the total physical memory usage */
std::uint64_t GetPhysicalMemoryUsage();

/* Get the total virtual memory usage */
std::uint64_t GetVirtualMemoryUsage();

} /* namespace MyGMapping */

#endif /* MY_GMAPPING_MEMORY_USAGE_HPP */
