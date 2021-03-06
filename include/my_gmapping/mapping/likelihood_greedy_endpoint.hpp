
/* likelihood_greedy_endpoint.hpp */

#ifndef MY_GMAPPING_LIKELIHOOD_GREEDY_ENDPOINT_HPP
#define MY_GMAPPING_LIKELIHOOD_GREEDY_ENDPOINT_HPP

#include "my_gmapping/mapping/likelihood_function.hpp"

namespace MyGMapping {
namespace Mapping {

class LikelihoodGreedyEndpoint final : public LikelihoodFunction
{
public:
    /* Constructor */
    LikelihoodGreedyEndpoint(const double mapResolution,
                             const double minUsableRange,
                             const double maxUsableRange,
                             const double hitAndMissedCellDist,
                             const double occupancyThreshold,
                             const double gaussianSigma,
                             const int kernelSize,
                             const double likelihoodScale);

    /* Destructor */
    ~LikelihoodGreedyEndpoint() = default;

    /* Calculate observation likelihood */
    double Likelihood(
        const GridMapInterface& gridMap,
        const Sensor::ScanDataPtr<double>& scanData,
        const RobotPose2D<double>& sensorPose) override;

private:
    /* Grid map resolution (in meters) */
    const double        mMapResolution;
    /* Minimum laser scan range considered for calculation */
    const double        mMinUsableRange;
    /* Maximum laser scan range considered for calculation */
    const double        mMaxUsableRange;
    /* Distance between the hit and the missed grid cell */
    const double        mHitAndMissedCellDist;
    /* Threshold probability to determine whether the grid cell is occupied */
    const double        mOccupancyThreshold;
    /* Threshold value to determine whether the grid cell is occupied */
    const std::uint16_t mOccupancyThresholdValue;
    /* Variance of the Gaussian distribution of the error */
    const double        mGaussianSigma;
    /* Size of the searching window (in the number of grid cells) */
    const int           mKernelSize;
    /* Scaling factor for the likelihood value */
    const double        mLikelihoodScale;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_LIKELIHOOD_GREEDY_ENDPOINT_HPP */
