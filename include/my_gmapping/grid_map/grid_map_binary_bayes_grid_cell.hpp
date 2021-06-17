
/* grid_map_binary_bayes_cell.hpp */

#ifndef MY_GMAPPING_GRID_MAP_GRID_MAP_BINARY_BAYES_CELL_HPP
#define MY_GMAPPING_GRID_MAP_GRID_MAP_BINARY_BAYES_CELL_HPP

#include "my_gmapping/grid_map/grid_map_cell.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <utility>

namespace MyGMapping {

template <typename T>
class BinaryBayesCell final : public Cell<T, T, double>
{
public:
    /* Base type */
    using BaseType = Cell<T, T, double>;

    /* Default constructor */
    BinaryBayesCell() : BaseType(),
                        mValue(Unknown) { }
    /* Destructor */
    ~BinaryBayesCell() = default;

    /* Copy constructor */
    BinaryBayesCell(const BinaryBayesCell&) = default;
    /* Copy assignment operator */
    BinaryBayesCell& operator=(const BinaryBayesCell&) = default;
    /* Move constructor */
    BinaryBayesCell(BinaryBayesCell&&) noexcept = default;
    /* Move assignment operator */
    BinaryBayesCell& operator=(BinaryBayesCell&&) noexcept = default;

    /* Cast operator */
    explicit operator T() const override { return this->mValue; }
    /* Get the occupancy probability value of the grid cell */
    T Value() const override { return this->mValue; }
    /* Get the internal representation of the grid cell */
    T RawValue() const override { return this->mValue; }

    /* Reset the grid cell state */
    void Reset() override;

    /* Set the value of the grid cell */
    void SetValue(const T newValue) override;
    /* Set the internal value of the grid cell */
    void SetRawValue(const T newRawValue) override;

    /* Update the value of the grid cell */
    void Update(const double probValue) override;

    /* Unknown occupancy probability value */
    using BaseType::Unknown;
    using BaseType::UnknownRaw;

    /* Smallest possible occupancy probability value */
    static constexpr T ProbabilityMin = static_cast<T>(1e-3);
    /* Largest possible occupancy probability value */
    static constexpr T ProbabilityMax = static_cast<T>(1.0 - 1e-3);

    /* Clamp a occupancy probability value */
    static T ClampValue(T probValue);
    /* Compute a odds from the occupancy probability value */
    static T ValueToOdds(T probValue);
    /* Compute an occupancy probability value from the odds */
    static T OddsToValue(T oddsValue);

private:
    /* Occupancy probability value */
    T mValue;
};

/* Reset the grid cell state */
template <typename T>
void BinaryBayesCell<T>::Reset()
{
    this->mValue = Unknown;
}

/* Set the value of the grid cell */
template <typename T>
void BinaryBayesCell<T>::SetValue(const T newValue)
{
    this->mValue = ClampValue(newValue);
}

/* Set the internal value of the grid cell */
template <typename T>
void BinaryBayesCell<T>::SetRawValue(const T newRawValue)
{
    this->mValue = ClampValue(newRawValue);
}

/* Update the value of the grid cell */
template <typename T>
void BinaryBayesCell<T>::Update(const double probValue)
{
    /* Assign a predefined probability value
     * if this grid cell is formerly not observed */
    if (this->mValue == Unknown) {
        /* Clamp the occupancy probability value for safety */
        this->mValue = ClampValue(probValue);
        return;
    }
    
    /* Otherwise, update the probability value using Binary Bayes Filter */
    const T oldOdds = ValueToOdds(this->mValue);
    const T valueOdds = ValueToOdds(probValue);
    const T newValue = OddsToValue(oldOdds * valueOdds);

    /* Clamp the occupancy probability value for safety */
    this->mValue = ClampValue(newValue);
}

/* Clamp a occupancy probability value */
template <typename T>
T BinaryBayesCell<T>::ClampValue(T probValue)
{
    /* Clamp a occupancy probability value */
    return std::clamp(probValue, ProbabilityMin, ProbabilityMax);
}

/* Compute a odds from the occupancy probability value */
template <typename T>
T BinaryBayesCell<T>::ValueToOdds(T probValue)
{
    /* Ensure that the probability value is not unknown */
    assert(probValue != Unknown);
    /* Clamp the occupancy probability value to avoid zero division */
    const T clampedValue = ClampValue(probValue);
    /* Compute the odds from the occupancy probability value */
    return clampedValue / (1.0 - clampedValue);
}

/* Compute an occupancy probability value from the odds */
template <typename T>
T BinaryBayesCell<T>::OddsToValue(T oddsValue)
{
    return ClampValue(oddsValue / (1.0 + oddsValue));
}

} /* namespace MyGMapping */

#endif /* MY_GMAPPING_GRID_MAP_GRID_MAP_BINARY_BAYES_CELL_HPP */
