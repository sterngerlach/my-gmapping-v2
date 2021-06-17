
/* grid_map_cell.hpp */

#ifndef MY_GMAPPING_GRID_MAP_GRID_MAP_CELL_HPP
#define MY_GMAPPING_GRID_MAP_GRID_MAP_CELL_HPP

#include <cstdint>

namespace MyGMapping {

template <typename T, typename U, typename V>
class Cell
{
public:
    /* Type definition */
    using ValueType = T;
    using StorageType = U;
    using ObservationType = V;

    /* Default constructor */
    Cell() = default;
    /* Destructor */
    virtual ~Cell() = default;

    /* Copy constructor */
    Cell(const Cell& other) = default;
    /* Copy assignment operator */
    Cell& operator=(const Cell& other) = default;
    /* Move constructor */
    Cell(Cell&& other) noexcept = default;
    /* Move assignment operator */
    Cell& operator=(Cell&& other) noexcept = default;

    /* Cast operator */
    virtual explicit operator T() const = 0;
    /* Get the value of the cell */
    virtual T Value() const = 0;
    /* Get the internal value of the grid cell */
    virtual U RawValue() const = 0;

    /* Reset the cell state */
    virtual void Reset() = 0;

    /* Set the value of the grid cell */
    virtual void SetValue(const T newValue) = 0;
    /* Set the internal value of the grid cell */
    virtual void SetRawValue(const U newRawValue) = 0;

    /* Update the value of the cell */
    virtual void Update(const V currentObservation) = 0;

    /* Unknown occupancy probability value */
    static constexpr T Unknown = static_cast<T>(0);
    /* Unknown occupancy probability value for internal representation */
    static constexpr U UnknownRaw = static_cast<U>(0);
};

} /* namespace MyGMapping */

#endif /* MY_GMAPPING_GRID_MAP_GRID_MAP_CELL_HPP */
