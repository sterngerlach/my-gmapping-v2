
/* scan_matcher_hill_climbing.hpp */

#ifndef MY_GMAPPING_MAPPING_SCAN_MATCHER_HILL_CLIMBING_HPP
#define MY_GMAPPING_MAPPING_SCAN_MATCHER_HILL_CLIMBING_HPP

#include "my_gmapping/mapping/likelihood_function.hpp"
#include "my_gmapping/mapping/scan_matcher.hpp"
#include "my_gmapping/metric/metric.hpp"

namespace MyGMapping {
namespace Mapping {

struct ScanMatcherHillClimbingMetrics
{
    /* Constructor */
    ScanMatcherHillClimbingMetrics(const std::string& scanMatcherName);
    /* Destructor */
    ~ScanMatcherHillClimbingMetrics() = default;

    /* Total processing time for the optimization */
    Metric::DistributionBase*         mOptimizationTime;
    /* Distance between the initial pose and the final pose */
    Metric::DistributionBase*         mDiffTranslation;
    /* Absolute difference between the initial angle and the final angle */
    Metric::DistributionBase*         mDiffRotation;
    /* Total number of the iterations */
    Metric::DistributionBase*         mNumOfIterations;
    /* Total number of the step size updates */
    Metric::DistributionBase*         mNumOfRefinements;
    /* Normalized likelihood value of the best particle */
    Metric::ValueSequenceBase<float>* mLikelihoodValue;
    /* Number of the scan points in the given scan */
    Metric::ValueSequenceBase<int>*   mNumOfScans;
};

class ScanMatcherHillClimbing final : public ScanMatcher
{
public:
    /* Type definitions */
    using LikelihoodFuncPtr = std::unique_ptr<LikelihoodFunction>;

    /* Constructor */
    ScanMatcherHillClimbing(
        const std::string& scanMatcherName,
        LikelihoodFuncPtr&& likelihoodFunc,
        const double linearDelta,
        const double angularDelta,
        const int maxIterations,
        const int numOfRefinements);

    /* Destructor */
    ~ScanMatcherHillClimbing() = default;

    /* Optimize particle pose by scan matching methods */
    void OptimizePose(
        const std::size_t numOfParticles,
        const std::vector<const GridMapType*>& particleMaps,
        const Sensor::ScanDataPtr<double>& scanData,
        const std::vector<RobotPose2D<double>>& initialPoses,
        std::vector<RobotPose2D<double>>& estimatedPoses,
        std::vector<double>& likelihoodValues) override;

private:
    /* Optimize particle pose by scan matching methods */
    void OptimizePoseCore(
        const GridMapInterfaceType& gridMap,
        const Sensor::ScanDataPtr<double>& scanData,
        const RobotPose2D<double>& initialPose,
        RobotPose2D<double>& estimatedPose,
        double& likelihoodValue);

private:
    /* Likelihood function to compute the observation likelihood */
    LikelihoodFuncPtr              mLikelihoodFunc;
    /* Initial step of the linear components (x and y) */
    const double                   mLinearDelta;
    /* Initial step of the angular component (theta) */
    const double                   mAngularDelta;
    /* Maximum number of the iterations */
    const int                      mMaxIterations;
    /* Maximum number of the step parameter updates */
    const int                      mNumOfRefinements;
    /* Metrics information */
    ScanMatcherHillClimbingMetrics mMetrics;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_MAPPING_SCAN_MATCHER_HILL_CLIMBING_HPP */
