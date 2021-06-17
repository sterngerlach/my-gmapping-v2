
/* grid_map_const_grid_cell.hpp */

#ifndef MY_GMAPPING_GRID_MAP_GRID_MAP_CONST_GRID_CELL_HPP
#define MY_GMAPPING_GRID_MAP_GRID_MAP_CONST_GRID_CELL_HPP

#include "my_gmapping/grid_map/grid_map_cell.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <limits>

namespace MyGMapping {

class ConstCell final : public Cell<double, std::uint16_t, double>
{
public:
    /* Base type */
    using BaseType = Cell<double, std::uint16_t, double>;

    /* Default constructor */
    ConstCell() : BaseType(),
                  mValue(UnknownRaw) { }
    /* Destructor */
    ~ConstCell() = default;

    /* Default copy constructor */
    ConstCell(const ConstCell&) = default;
    /* Default copy assignment operator */
    ConstCell& operator=(const ConstCell&) = default;
    /* Default move constructor */
    ConstCell(ConstCell&&) noexcept = default;
    /* Default move assignment operator */
    ConstCell& operator=(ConstCell&&) noexcept = default;

    /* Cast operator */
    explicit operator double() const override
    { return RawToValue(this->mValue); }
    /* Get the value of the cell */
    double Value() const override
    { return RawToValue(this->mValue); }
    /* Get the internal representation of the cell */
    std::uint16_t RawValue() const override { return this->mValue; }

    /* Reset the cell state */
    void Reset() override { this->mValue = UnknownRaw; }

    /* Set the value of the grid cell */
    void SetValue(const double newValue) override
    { this->mValue = ValueToRaw(newValue); }
    /* Set the internal value of the grid cell */
    void SetRawValue(const std::uint16_t newRawValue) override
    { this->mValue = newRawValue; }

    /* Update the value of the cell */
    void Update(const double newValue) override
    { this->mValue = ValueToRaw(newValue); }

    /* Smallest possible occupancy probability value */
    static constexpr std::uint16_t ValueMin = 1;
    /* Largest possible occupancy probability value */
    static constexpr std::uint16_t ValueMax =
        std::numeric_limits<std::uint16_t>::max();
    /* Occupancy probability value range */
    static constexpr std::uint16_t ValueRange = ValueMax - ValueMin;

    /* Convert the internal representation to the probability value */
    static double RawToValue(const std::uint16_t discretizedValue);
    /* Convert the probability value to the internal representation */
    static std::uint16_t ValueToRaw(const double probValue);

private:
    /* Occupancy probability value */
    std::uint16_t mValue;
};

} /* namespace MyGMapping */

#endif /* MY_GMAPPING_GRID_MAP_GRID_MAP_CONST_GRID_CELL_HPP */
