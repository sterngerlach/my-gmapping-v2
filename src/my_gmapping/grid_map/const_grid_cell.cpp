
/* const_grid_cell.cpp */

#include "my_gmapping/grid_map/grid_map_const_grid_cell.hpp"

namespace MyGMapping {

/* Convert the internal representation to the probability value */
double ConstCell::RawToValue(const std::uint16_t discretizedValue)
{
    if (discretizedValue == UnknownRaw)
        return Unknown;

    return static_cast<double>(discretizedValue - ValueMin) /
           static_cast<double>(ValueRange);
}

/* Convert the probability value to the internal representation */
std::uint16_t ConstCell::ValueToRaw(const double probValue)
{
    /* Ensure that the probability value is valid */
    assert(probValue >= 0.0 && probValue <= 1.0);

    if (probValue == Unknown)
        return UnknownRaw;
    if (probValue < 0.0)
        return ValueMin;
    if (probValue > 1.0)
        return ValueMax;

    return static_cast<std::uint16_t>(probValue * ValueRange) + ValueMin;
}

} /* namespace MyGMapping */
