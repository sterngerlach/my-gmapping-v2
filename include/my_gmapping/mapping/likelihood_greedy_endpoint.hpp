
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
                             const double likelihoodScale,
                             const bool useLocalMap,
                             const int localMapSize);

    /* Destructor */
    ~LikelihoodGreedyEndpoint() = default;

    /* Calculate observation likelihood */
    double Likelihood(
        const GridMapInterfaceType& gridMap,
        const Sensor::ScanDataPtr<double>& scanData,
        const RobotPose2D<double>& sensorPose) override;

private:
    /* Setup the lookup table for the likelihood value */
    void SetupLookupTable();

private:
    /* Grid map resolution (in meters) */
    const double              mMapResolution;
    /* Minimum laser scan range considered for calculation */
    const double              mMinUsableRange;
    /* Maximum laser scan range considered for calculation */
    const double              mMaxUsableRange;
    /* Distance between the hit and the missed grid cell */
    const double              mHitAndMissedCellDist;
    /* Occupancy probability threshold for being obstructed */
    const double              mOccupancyThreshold;
    /* Variance of the Gaussian distribution of the error */
    const double              mGaussianSigma;
    /* Size of the searching window (in the number of grid cells) */
    const int                 mKernelSize;
    /* Scaling factor for the likelihood value */
    const double              mLikelihoodScale;
    /* Lookup table for the likelihood value */
    std::unique_ptr<double[]> mLikelihoodTable;
    /* Minimum likelihood value */
    double                    mDefaultLikelihood;

    /* TODO: Create a cropped version of a GridMap class and get rid of this */
    /* Flag to determine whether the cropped grid map is used */
    const bool                mUseLocalMap;
    /* Size of the cropped grid map (in the number of grid cells) */
    const int                 mLocalMapSize;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_LIKELIHOOD_GREEDY_ENDPOINT_HPP */
