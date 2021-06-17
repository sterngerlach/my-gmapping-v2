
/* grid_map_types.hpp */

#ifndef MY_GMAPPING_MAPPING_GRID_MAP_TYPES_HPP
#define MY_GMAPPING_MAPPING_GRID_MAP_TYPES_HPP

#include "my_gmapping/grid_map_new/grid_binary_bayes.hpp"
#include "my_gmapping/grid_map_new/grid_constant.hpp"
#include "my_gmapping/grid_map_new/grid_map.hpp"

#include <type_traits>

namespace MyGMapping {
namespace Mapping {

/* Make sure that grid maps and precomputed grid maps store and represent
 * occupancy probability values using the same data type */
static_assert(std::is_same<GridMapNew::GridBinaryBayes::ProbabilityType,
                           GridMapNew::GridConstant::ProbabilityType>::value,
              "GridBinaryBayes and GridConstant should represent "
              "occupancy probability values using the same data type");
static_assert(std::is_same<GridMapNew::GridBinaryBayes::ValueType,
                           GridMapNew::GridConstant::ValueType>::value,
              "GridBinaryBayes and GridConstant should store "
              "occupancy probability values using the same data type");

/* Type definitions */
using GridMap = GridMapNew::GridMap<GridMapNew::GridBinaryBayes>;
using ConstMap = GridMapNew::GridMap<GridMapNew::GridConstant>;
using GridMapInterface = GridMapNew::GridMapInterface<
    GridMapNew::GridBinaryBayes::ProbabilityType,
    GridMapNew::GridBinaryBayes::ValueType>;

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_MAPPING_GRID_MAP_TYPES_HPP */
