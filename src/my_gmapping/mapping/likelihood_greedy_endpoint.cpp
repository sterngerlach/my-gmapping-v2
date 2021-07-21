
/* likelihood_greedy_endpoint.cpp */

#include "my_gmapping/mapping/likelihood_greedy_endpoint.hpp"

namespace MyGMapping {
namespace Mapping {

/* Constructor */
LikelihoodGreedyEndpoint::LikelihoodGreedyEndpoint(
    const double mapResolution,
    const double minUsableRange,
    const double maxUsableRange,
    const double hitAndMissedCellDist,
    const double occupancyThreshold,
    const double gaussianSigma,
    const int kernelSize,
    const double likelihoodScale) :
    LikelihoodFunction(),
    mMapResolution(mapResolution),
    mMinUsableRange(minUsableRange),
    mMaxUsableRange(maxUsableRange),
    mHitAndMissedCellDist(hitAndMissedCellDist),
    mOccupancyThreshold(occupancyThreshold),
    mOccupancyThresholdValue(
        GridMap::GridType::ProbabilityToValue(occupancyThreshold)),
    mGaussianSigma(gaussianSigma),
    mKernelSize(kernelSize),
    mLikelihoodScale(likelihoodScale)
{
}

/* Calculate observation likelihood based on squared distance
 * between scan point and matched cell position */
double LikelihoodGreedyEndpoint::Likelihood(
    const GridMapInterface& gridMap,
    const Sensor::ScanDataPtr<double>& scanData,
    const RobotPose2D<double>& sensorPose)
{
    /* Sum of the squared distances between hit points and
     * their corresponding grid cells */
    int distSquaredSum = 0;
    /* Compute the maximum squared distance */
    const int distSquaredMax =
        2 * (this->mKernelSize + 1) * (this->mKernelSize + 1);

    const double minRange = std::max(
        this->mMinUsableRange, scanData->MinRange());
    const double maxRange = std::min(
        this->mMaxUsableRange, scanData->MaxRange());

    for (std::size_t i = 0; i < scanData->NumOfScans(); ++i) {
        const double range = scanData->RangeAt(i);

        if (range >= maxRange || range <= minRange)
            continue;

        /* Compute the grid cell indices corresponding to the scan point and
         * the missed point */
        Point2D<double> hitPoint;
        Point2D<double> missedPoint;
        scanData->HitAndMissedPoint(
            sensorPose, i, this->mHitAndMissedCellDist,
            hitPoint, missedPoint);

        const Point2D<int> hitIdx =
            gridMap.PositionToIndex(hitPoint.mX, hitPoint.mY);
        const Point2D<int> missedIdx =
            gridMap.PositionToIndex(missedPoint.mX, missedPoint.mY);

        /* Search the best cell index from the window */
        const std::uint16_t unknownValue = gridMap.UnknownValue();
        /* Initialize the minimum squared distance with the default value */
        int distSquaredMin = distSquaredMax;

        for (int ky = -this->mKernelSize; ky <= this->mKernelSize; ++ky) {
            for (int kx = -this->mKernelSize; kx <= this->mKernelSize; ++kx) {
                const std::uint16_t hitValue = gridMap.ValueOr(
                    hitIdx.mY + ky, hitIdx.mX + kx, unknownValue);
                const std::uint16_t missedValue = gridMap.ValueOr(
                    missedIdx.mY + ky, missedIdx.mX + kx, unknownValue);

                /* Skip if the grid cell is not observed yet */
                if (hitValue == unknownValue || missedValue == unknownValue)
                    continue;
                /* Skip if the value of the hit grid cell is less than the
                 * occupancy threshold or the value of the missed one is
                 * greater than the occupancy threshold */
                if (hitValue < this->mOccupancyThresholdValue ||
                    missedValue > this->mOccupancyThresholdValue)
                    continue;

                /* Compute the squared distance between the hit point and
                 * its corresponding grid cell */
                const int distSquared = kx * kx + ky * ky;
                /* Update the minimum squared distance */
                distSquaredMin = std::min(distSquaredMin, distSquared);
            }
        }

        /* Update the sum of squared distances */
        distSquaredSum += distSquaredMin;
    }

    /* Compute the matching score from the sum of squared distances
     * which represents the log-likelihood of the observation */
    double likelihoodSum =
        (-this->mMapResolution * this->mMapResolution * distSquaredSum) /
        (0.5 * this->mGaussianSigma);
    /* Scale the matching score to prevent underflow in weight computation */
    likelihoodSum *= this->mLikelihoodScale;

    return likelihoodSum;
}

} /* namespace Mapping */
} /* namespace MyGMapping */
