
/* grid_map_counting_grid_cell.hpp */

#ifndef MY_GMAPPING_GRID_MAP_GRID_MAP_COUNTING_GRID_CELL_HPP
#define MY_GMAPPING_GRID_MAP_GRID_MAP_COUNTING_GRID_CELL_HPP

#include "my_gmapping/grid_map/grid_map_cell.hpp"
#include "my_gmapping/util.hpp"

#include <cstdint>
#include <utility>

namespace MyGMapping {

template <typename T>
class CountingCell final : public Cell<T, T, bool>
{
public:
    /* Base type */
    using BaseType = Cell<T, T, bool>;

    /* Default constructor */
    CountingCell() : BaseType(),
                     mValue(Unknown),
                     mHit(0),
                     mMiss(0) { }
    /* Destructor */
    ~CountingCell() = default;

    /* Copy constructor */
    CountingCell(const CountingCell& other) = default;
    /* Copy assignment operator */
    CountingCell& operator=(const CountingCell& other) = default;
    /* Move constructor */
    CountingCell(CountingCell&& other) noexcept = default;
    /* Move assignment operator */
    CountingCell& operator=(CountingCell&& other) noexcept = default;

    /* Cast operator */
    explicit operator T() const override { return this->mValue; }
    /* Get the value of the cell */
    T Value() const override { return this->mValue; }
    /* Get the internal representation of the grid cell */
    T RawValue() const override { return this->mValue; }

    /* Reset the cell state */
    void Reset() override;

    /* Set the value of the grid cell (Not implemented) */
    void SetValue(const T newValue) override
    { XAssert(false, "SetValue() is not implemented"); }
    /* Set the internal value of the grid cell (Not implemented) */
    void SetRawValue(const T newRawValue) override
    { XAssert(false, "SetRawValue() is not implemented"); }

    /* Update the value of the cell */
    void Update(const bool hitOrMiss) override;

    /* Unknown occupancy probability value */
    using BaseType::Unknown;
    using BaseType::UnknownRaw;

    /* Smallest possible occupancy probability value,
     * which is slightly larger than the unknown value */
    static constexpr T Epsilon = static_cast<T>(1e-3);

private:
    T             mValue;
    std::uint32_t mHit;
    std::uint32_t mMiss;
};

/* Reset the cell state */
template <typename T>
void CountingCell<T>::Reset()
{
    this->mValue = Unknown;
    this->mHit = 0;
    this->mMiss = 0;
}

/* Update the value of the cell */
template <typename T>
void CountingCell<T>::Update(const bool hitOrMiss)
{
    /* Update the counter */
    if (hitOrMiss)
        ++this->mHit;
    else
        ++this->mMiss;

    /* Update the occupancy probability value */
    this->mValue = static_cast<T>(this->mHit) /
                   static_cast<T>(this->mHit + this->mMiss);

    /* If the number of the hits is zero, then add the small value
     * to distinguish from the unknown state */
    this->mValue = !this->mHit ? Epsilon : this->mValue;
}

} /* namespace MyGMapping */

#endif /* MY_GMAPPING_GRID_MAP_GRID_MAP_COUNTING_GRID_CELL_HPP */
