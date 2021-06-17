
/* weight_normalizer.hpp */

#ifndef MY_GMAPPING_WEIGHT_NORMALIZER_HPP
#define MY_GMAPPING_WEIGHT_NORMALIZER_HPP

#include <cmath>
#include <vector>

namespace MyGMapping {
namespace Mapping {

/*
 * WeightNormalizationType enum specifies the method of weight normalization
 */
enum class WeightNormalizationType
{
    Basic,
    ExponentialWeight,
};

/*
 * WeightNormalizer class is for normalizing particle weights
 */
class WeightNormalizer final
{
public:
    /* Constructor */
    WeightNormalizer() = default;
    /* Destructor */
    ~WeightNormalizer() = default;

    /* Normalize weights */
    void NormalizeWeights(
        const WeightNormalizationType normalizationType,
        std::vector<double>& particleWeights) const;

private:
    /* Normalize weights (basic version) */
    void NormalizeWeightsBasic(
        std::vector<double>& particleWeights) const;

    /* Normalize weights (exponentially weighted version) */
    void NormalizeWeightsExponential(
        std::vector<double>& particleWeights) const;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_WEIGHT_NORMALIZER_HPP */
