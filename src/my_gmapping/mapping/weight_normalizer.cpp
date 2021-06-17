
/* weight_normalizer.cpp */

#include "my_gmapping/mapping/weight_normalizer.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>

namespace MyGMapping {
namespace Mapping {

/* Normalize weights */
void WeightNormalizer::NormalizeWeights(
    const WeightNormalizationType normalizationType,
    std::vector<double>& particleWeights) const
{
    switch (normalizationType) {
        case WeightNormalizationType::Basic:
            return this->NormalizeWeightsBasic(particleWeights);
        case WeightNormalizationType::ExponentialWeight:
            return this->NormalizeWeightsExponential(particleWeights);
        default:
            /* Unknown particle weight normalization method is specified */
            assert(false);
            break;
    }

    return;
}

/* Normalize weights (basic version) */
void WeightNormalizer::NormalizeWeightsBasic(
    std::vector<double>& particleWeights) const
{
    /* Compute the sum of the particle weights */
    const double weightSum = std::accumulate(
        particleWeights.cbegin(), particleWeights.cend(), 0.0);

    /* Normalize the particle weights */
    std::transform(
        particleWeights.begin(), particleWeights.end(),
        particleWeights.begin(),
        [weightSum](double weight) { return weight / weightSum; });
}

/* Normalize weights (exponentially weighted version) */
void WeightNormalizer::NormalizeWeightsExponential(
    std::vector<double>& particleWeights) const
{
    /* Subtract by the max weight to avoid overflows */
    const auto maxIter = std::max_element(
        particleWeights.cbegin(), particleWeights.cend());
    const double maxWeight = *maxIter;

    std::transform(
        particleWeights.begin(), particleWeights.end(),
        particleWeights.begin(),
        [maxWeight](double weight) { return std::exp(weight - maxWeight); });

    /* Compute the sum of the particle weights */
    const double weightSum = std::accumulate(
        particleWeights.cbegin(), particleWeights.cend(), 0.0);

    /* Normalize the particle weights */
    std::transform(
        particleWeights.begin(), particleWeights.end(),
        particleWeights.begin(),
        [weightSum](double weight) { return weight / weightSum; });
}

} /* namespace Mapping */
} /* namespace MyGMapping */
