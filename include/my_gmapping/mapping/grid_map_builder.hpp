
/* grid_map_builder.hpp */

#ifndef MY_GMAPPING_GRID_MAP_BUILDER_HPP
#define MY_GMAPPING_GRID_MAP_BUILDER_HPP

#include <deque>
#include <memory>
#include <random>
#include <vector>

#include "my_gmapping/point.hpp"
#include "my_gmapping/pose.hpp"
#include "my_gmapping/mapping/covariance_estimator.hpp"
#include "my_gmapping/mapping/grid_map_types.hpp"
#include "my_gmapping/mapping/likelihood_function.hpp"
#include "my_gmapping/mapping/map_builder.hpp"
#include "my_gmapping/mapping/motion_model.hpp"
#include "my_gmapping/mapping/particle.hpp"
#include "my_gmapping/mapping/scan_interpolator.hpp"
#include "my_gmapping/mapping/scan_matcher.hpp"
#include "my_gmapping/mapping/trajectory_node.hpp"
#include "my_gmapping/mapping/weight_normalizer.hpp"
#include "my_gmapping/metric/metric.hpp"
#include "my_gmapping/sensor/sensor_data.hpp"

namespace MyGMapping {
namespace Mapping {

struct GridMapBuilderMetrics
{
    /* Constructor */
    GridMapBuilderMetrics();
    /* Destructor */
    ~GridMapBuilderMetrics() = default;

    /* Total number of the input scan data */
    Metric::CounterBase*                      mInputScanDataCount;
    /* Total number of the processed scan data */
    Metric::CounterBase*                      mProcessCount;
    /* Total processing time of the SLAM frontend */
    Metric::ValueSequenceBase<int>*           mProcessTime;
    /* Total processing time for the scan data */
    Metric::ValueSequenceBase<int>*           mProcessScanTime;
    /* Total processing time for the particle sampling */
    Metric::ValueSequenceBase<int>*           mSamplingTime;
    /* Total processing time for setting up the scan data */
    Metric::ValueSequenceBase<int>*           mScanDataSetupTime;
    /* Total processing time for the scan matching */
    Metric::ValueSequenceBase<int>*           mScanMatchingTime;
    /* Total processing time for the final scan matching */
    Metric::ValueSequenceBase<int>*           mFinalScanMatchingTime;
    /* Total processing time for updating the particle weights */
    Metric::ValueSequenceBase<int>*           mWeightUpdateTime;
    /* Total processing time for updating the grid map */
    Metric::ValueSequenceBase<int>*           mMapUpdateTime;
    /* Total processing time for updating the latest grid map */
    Metric::ValueSequenceBase<int>*           mLatestMapUpdateTime;
    /* Total processing time for the particle resampling */
    Metric::ValueSequenceBase<int>*           mResamplingTime;
    /* Accumulated travel distance between the processed scans */
    Metric::ValueSequenceBase<float>*         mIntervalTravelDist;
    /* Difference of the robot pose angle between the processed scans */
    Metric::ValueSequenceBase<float>*         mIntervalAngle;
    /* Time between the processed scans */
    Metric::ValueSequenceBase<float>*         mIntervalTime;
    /* Number of the scan points for each scan data */
    Metric::ValueSequenceBase<int>*           mNumOfScans;
    /* Frame number of the processed scan data */
    Metric::ValueSequenceBase<int>*           mProcessFrame;
    /* Effective sample size */
    Metric::ValueSequenceBase<float>*         mEffectiveSampleSize;
    /* Total memory consumption for the particle grid maps */
    Metric::ValueSequenceBase<std::uint64_t>* mGridMapMemoryUsage;
    /* Total memory consumption for the latest grid maps */
    Metric::ValueSequenceBase<std::uint64_t>* mLatestMapMemoryUsage;
    /* Total memory consumption for the trajectory nodes */
    Metric::ValueSequenceBase<std::uint64_t>* mTrajectoryMemoryUsage;
};

class GridMapBuilder
{
public:
    /* Type definition */
    using ScanDataDeque = std::deque<Sensor::ScanDataPtr<double>>;

    /* Constructor */
    GridMapBuilder(
        std::unique_ptr<MotionModel>&& motionModel,
        std::unique_ptr<LikelihoodFunction>&& likelihoodFunc,
        std::unique_ptr<ScanMatcher>&& scanMatcher,
        std::unique_ptr<ScanMatcher>&& finalScanMatcher,
        std::unique_ptr<ScanInterpolator>&& scanInterpolator,
        std::unique_ptr<CovarianceEstimator>&& covarianceEstimator,
        MapBuilder&& mapBuilder,
        const WeightNormalizationType weightNormalizationType,
        const int numOfParticles,
        const bool useLatestMap,
        const std::size_t numOfScansForLatestMap,
        const RobotPose2D<double>& initialPose,
        const double mapResolution,
        const double updateThresholdTravelDist,
        const double updateThresholdAngle,
        const double updateThresholdTime,
        const double resampleThreshold,
        const double degenerationThreshold);
    /* Destructor */
    ~GridMapBuilder() = default;

    /* Copy constructor */
    GridMapBuilder(const GridMapBuilder& other) = delete;
    /* Copy assignment operator */
    GridMapBuilder& operator=(const GridMapBuilder& other) = delete;
    /* Move constructor */
    GridMapBuilder(GridMapBuilder&& other) = delete;
    /* Move assignment operator */
    GridMapBuilder& operator=(GridMapBuilder&& other) = delete;

    /* Process scan data */
    bool ProcessScan(const Sensor::ScanDataPtr<double>& rawScanData,
                     const RobotPose2D<double>& odomPose);

    /* Get the process counter */
    inline int ProcessCounter() const { return this->mProcessCount; }

    /* Get the number of particles */
    inline int NumOfParticles() const
    { return static_cast<int>(this->mParticles.size()); }

    /* Get the map of the specified particle */
    inline const Particle& ParticleAt(std::size_t particleIdx) const
    { return this->mParticles[particleIdx]; }

    /* Get the index of the best particle */
    std::size_t BestParticleIndex() const;

    /* Get the estimated trajectory of the specified particle */
    std::vector<TimeStampedPose> ParticleTrajectoryWithTimeStamp(
        std::size_t particleIdx) const;

    /* Get the estimated trajectory of the specified particle */
    std::vector<RobotPose2D<double>> ParticleTrajectory(
        std::size_t particleIdx) const;

private:
    /* Execute scan matching and update particle weights */
    void ExecuteScanMatching(const Sensor::ScanDataPtr<double>& scanData);

    /* Update the deque that stores the latest scan data */
    void UpdateLatestScans(const Sensor::ScanDataPtr<double>& scanData);

    /* Update the grid maps for all particles */
    void UpdateGridMaps(const Sensor::ScanDataPtr<double>& scanData);
    /* Update the latest maps for all particles */
    void UpdateLatestMaps(const ScanDataDeque& latestScanData);

    /* Resample particles according to their weights if necessary */
    void ResampleParticles();

    /* Calculate the effective sample size of the particles */
    double CalculateEffectiveSampleSize() const;

    /* Check the degeneration */
    bool CheckDegeneration(const Eigen::Matrix3d& poseCovarianceMat) const;

    /* Inspect the memory usage for the particle grid maps */
    std::uint64_t InspectGridMapMemoryUsage() const;
    /* Inspect the memory usage for the latest grid maps */
    std::uint64_t InspectLatestMapMemoryUsage() const;
    /* Inspect the memory usage for the trajectory nodes */
    std::uint64_t InspectTrajectoryMemoryUsage() const;

private:
    /* Total number of the processed input data */
    int                     mProcessCount;
    /* Resolution of the grid map */
    const double            mResolution;
    /* Collection of the particles */
    std::vector<Particle>   mParticles;
    /* Motion model to sample the particle poses using the odometry */
    MotionModelPtr          mMotionModel;
    /* Likelihood function to compute the observation likelihood */
    LikelihoodFunctionPtr   mLikelihoodFunc;
    /* Scan matcher to refine the particle poses */
    ScanMatcherPtr          mScanMatcher;
    /* Final scan matcher that executes the scan matching at sub-pixel
     * accuracy and refines the result returned from the first scan matcher */
    ScanMatcherPtr          mFinalScanMatcher;
    /* Scan interpolator to correct the raw scan data */
    ScanInterpolatorPtr     mScanInterpolator;
    /* Covariance estimator to compute the covariance of particle pose */
    CovarianceEstimatorPtr  mCovarianceEstimator;
    /* Map builder to update the grid map using the scan data */
    MapBuilder              mMapBuilder;
    /* Weight normalizer to normalize the particle importance weights */
    WeightNormalizer        mWeightNormalizer;
    /* Method of the weight normalization */
    WeightNormalizationType mWeightNormalizationType;
    /* Random engine */
    std::mt19937            mRandEngine;
    /* Flag to determine whether the latest map is used for scan matching */
    bool                    mUseLatestMap;
    /* The number of scans used to construct the latest map */
    std::size_t             mNumOfScansForLatestMap;
    /* Deque of the latest scan data for constructing the latest map */
    ScanDataDeque           mLatestScanData;

    /* Last odometry pose */
    RobotPose2D<double>     mLastOdomPose;
    /* Accumulated travel distance since the last map update */
    double                  mAccumulatedTravelDist;
    /* Accumulated angle since the last map update */
    double                  mAccumulatedAngle;
    /* Odometry pose at the last map update */
    RobotPose2D<double>     mLastMapUpdateOdomPose;
    /* Time of the last map update */
    double                  mLastMapUpdateTime;
    /* Map update threshold for accumulated travel distance */
    const double            mUpdateThresholdTravelDist;
    /* Map update threshold for accumulated angle */
    const double            mUpdateThresholdAngle;
    /* Map update threshold for the elapsed time since the last map update */
    const double            mUpdateThresholdTime;
    /* Effective sample size threshold for the particle sampling */
    const double            mResampleThreshold;
    /* Eigenvalues threshold for checking degeneration */
    const double            mDegenerationThreshold;
    /* Metrics information */
    GridMapBuilderMetrics   mMetrics;
};

/*
 * Utility function declarations
 */

/* Compute the maximum of a 'winSize' pixel wide row at each pixel */
void SlidingWindowMaxRow(const GridMap& gridMap,
                         const int winSize,
                         const Point2D<int>& idxMin,
                         ConstMap& intermediateMap);

/* Compute the maximum of a 'winSize' pixel wide column at each pixel */
void SlidingWindowMaxCol(const ConstMap& intermediateMap,
                         const int winSize,
                         ConstMap& precompMap);

/* Precompute grid map for efficiency */
void PrecomputeGridMap(const GridMap& gridMap,
                       const int winSize,
                       const Point2D<int>& idxMin,
                       ConstMap& intermediateMap,
                       ConstMap& precompMap);

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_GRID_MAP_BUILDER_HPP */
