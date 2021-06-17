
/* scan_matcher_correlative.hpp */

#ifndef MY_GMAPPING_MAPPING_SCAN_MATCHER_CORRELATIVE_HPP
#define MY_GMAPPING_MAPPING_SCAN_MATCHER_CORRELATIVE_HPP

#include "my_gmapping/mapping/scan_matcher.hpp"

#include "my_gmapping/grid_map/grid_map.hpp"
#include "my_gmapping/mapping/grid_map_types.hpp"
#include "my_gmapping/mapping/likelihood_function.hpp"
#include "my_gmapping/metric/metric.hpp"

namespace MyGMapping {
namespace Mapping {

struct ScanMatcherCorrelativeMetrics
{
    /* Constructor */
    ScanMatcherCorrelativeMetrics(const std::string& scanMatcherName);
    /* Destructor */
    ~ScanMatcherCorrelativeMetrics() = default;

    /* Total processing time for setting up the input */
    Metric::DistributionBase*         mInputSetupTime;
    /* Total processing time for the optimization */
    Metric::DistributionBase*         mOptimizationTime;
    /* Distance between the initial pose and the final pose */
    Metric::DistributionBase*         mDiffTranslation;
    /* Absolute difference between the initial angle and the final angle */
    Metric::DistributionBase*         mDiffRotation;
    /* Size of the search window along the X-axis */
    Metric::DistributionBase*         mWinSizeX;
    /* Size of the search window along the Y-axis */
    Metric::DistributionBase*         mWinSizeY;
    /* Size of the search window along the Theta-axis */
    Metric::DistributionBase*         mWinSizeTheta;
    /* Step size along the X-axis */
    Metric::DistributionBase*         mStepSizeX;
    /* Step size along the Y-axis */
    Metric::DistributionBase*         mStepSizeY;
    /* Step size along the Theta-axis */
    Metric::DistributionBase*         mStepSizeTheta;
    /* Total number of the skipped nodes */
    Metric::CounterBase*              mNumOfIgnoredNodes;
    /* Total number of the processed nodes */
    Metric::CounterBase*              mNumOfProcessedNodes;
    /* Normalized score value of the best particle */
    Metric::ValueSequenceBase<float>* mScoreValue;
    /* Normalized likelihood value of the best particle */
    Metric::ValueSequenceBase<float>* mLikelihoodValue;
    /* Number of the scan points in the given scan */
    Metric::ValueSequenceBase<int>*   mNumOfScans;
};

class ScanMatcherCorrelative final : public ScanMatcher
{
public:
    /* Type definitions */
    using LikelihoodFuncPtr = std::unique_ptr<LikelihoodFunction>;

    /* Constructor */
    ScanMatcherCorrelative(const std::string& scanMatcherName,
                           LikelihoodFuncPtr&& likelihoodFunc,
                           const bool useCroppedMap,
                           const int croppedMapSizeX,
                           const int croppedMapSizeY,
                           const int lowResolution,
                           const double rangeX,
                           const double rangeY,
                           const double rangeTheta,
                           const double scanRangeMax);
    /* Destructor */
    ~ScanMatcherCorrelative() = default;

    /* Optimize the particle poses based on the correlative scan matching */
    ScanMatchingResultVector OptimizePose(
        const ScanMatchingQueryVector& queries,
        const Sensor::ScanDataPtr<double>& scanData) override;

private:
    /* Optimize the particle pose based on the correlative scan matching */
    ScanMatchingResult OptimizePoseCore(
        const GridMapInterface& gridMap,
        const BoundingBox<int>& boundingBox,
        const GridMapInterface& coarseGridMap,
        const Sensor::ScanDataPtr<double>& scanData,
        const RobotPose2D<double>& initialPose,
        const double normalizedScoreThreshold);

    /* Compute the search step */
    void ComputeSearchStep(const GridMapInterface& gridMap,
                           const Sensor::ScanDataPtr<double>& scanData,
                           double& stepX,
                           double& stepY,
                           double& stepTheta) const;

    /* Compute the grid cell indices for scan points */
    void ComputeScanIndices(
        const GridMapInterface& gridMap,
        const RobotPose2D<double>& sensorPose,
        const Sensor::ScanDataPtr<double>& scanData,
        std::vector<Point2D<int>>& scanIndices) const;

    /* Compute the normalized scan matching score based on the
     * already projected scan points (indices) and index offsets */
    double ComputeScore(
        const GridMapInterface& gridMap,
        const BoundingBox<int>& boundingBox,
        const std::vector<Point2D<int>>& scanIndices,
        const Point2D<int>& scanIdxOffset,
        const int offsetX,
        const int offsetY) const;

    /* Evaluate the matching score using high-resolution grid map */
    void EvaluateHighResolutionMap(
        const GridMapInterface& gridMap,
        const BoundingBox<int>& boundingBox,
        const std::vector<Point2D<int>>& scanIndices,
        const Point2D<int>& scanIdxOffset,
        const int offsetX,
        const int offsetY,
        const int offsetTheta,
        int& maxWinX,
        int& maxWinY,
        int& maxWinTheta,
        double& maxScore) const;

private:
    /* Likelihood function to compute the observation likelihood */
    LikelihoodFuncPtr             mLikelihoodFunc;
    /* Flag to determine whether the grid map is cropped */
    const bool                    mUseCroppedMap;
    /* Width of the cropped grid map (in the number of grid cells) */
    const int                     mCroppedMapSizeX;
    /* Height of the cropped grid map (in the number of grid cells) */
    const int                     mCroppedMapSizeY;
    /* Resolution for low resolution map (in the number of grid cells) */
    const int                     mLowResolution;
    /* Width of the searching window (in meters) */
    const double                  mRangeX;
    /* Height of the searching window (in meters) */
    const double                  mRangeY;
    /* Angular range of the searching window (in radians) */
    const double                  mRangeTheta;
    /* Maximum laser scan range considered for scan matching */
    const double                  mScanRangeMax;
    /* Metrics information */
    ScanMatcherCorrelativeMetrics mMetrics;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_MAPPING_SCAN_MATCHER_CORRELATIVE_HPP */
