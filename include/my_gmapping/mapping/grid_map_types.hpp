
/* grid_map_types.hpp */

#ifndef MY_GMAPPING_MAPPING_GRID_MAP_TYPES_HPP
#define MY_GMAPPING_MAPPING_GRID_MAP_TYPES_HPP

#include "my_gmapping/grid_map/grid_map.hpp"
#include "my_gmapping/grid_map/grid_map_discrete_grid_cell.hpp"
#include "my_gmapping/grid_map/grid_map_const_grid_cell.hpp"

#include <type_traits>

namespace MyGMapping {
namespace Mapping {

/* Make sure that grid maps and precomputed grid maps store and represent
 * occupancy probability values using the same data type */
static_assert(std::is_same<DiscreteCell::ValueType,
                           ConstCell::ValueType>::value,
              "DiscreteCell and ConstCell should represent "
              "occupancy probability values using the same data type");
static_assert(std::is_same<DiscreteCell::StorageType,
                           ConstCell::StorageType>::value,
              "DiscreteCell and ConstCell should store "
              "occupancy probability values using the same data type");

/* Type declarations for the grid map and the precomputed grid map */
using GridMapType = GridMap<DiscreteCell>;
using ConstMapType = GridMap<ConstCell>;
using GridMapInterfaceType = GridMapBase<DiscreteCell::ValueType,
                                         DiscreteCell::StorageType>;

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_MAPPING_GRID_MAP_TYPES_HPP */
