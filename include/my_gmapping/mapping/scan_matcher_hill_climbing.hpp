
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
    /* `Times` struct holds the processing times for each particle */
    struct Times
    {
        int mOptimizationTime;
    };

    /* `Parameters` struct holds the parameter settings for each particle */
    struct Parameters
    {
        /* Metric values related to the algorithmic parameters */
        float mDiffTranslation;
        float mDiffRotation;
        /* Metric values related to the algorithmic performances */
        int   mNumOfIterations;
        int   mNumOfRefinements;
        /* Metric values related to the outputs */
        float mLikelihoodValue;
        int   mNumOfScans;
    };

    /* Constructor */
    ScanMatcherHillClimbingMetrics(const std::string& scanMatcherName);
    /* Destructor */
    ~ScanMatcherHillClimbingMetrics() = default;

    /* Resize the buffer to store the processing times */
    void Resize(const std::size_t numOfParticles);
    /* Set the processing times for each particle */
    void SetTimes(const std::size_t idx, const Times& times);
    /* Set the parameter settings for each particle */
    void SetParameters(const std::size_t idx, const Parameters& params);
    /* Collect the particle-wise metrics and update the overall metrics */
    void Update(const std::size_t bestIdx);

    /* Processing times for all particles */
    std::vector<Times>                mTimes;
    /* Parameter settings for all particles */
    std::vector<Parameters>           mParams;

    /* Total processing time for the optimization */
    Metric::ValueSequenceBase<int>*   mOptimizationTime;
    /* Distance between the initial pose and the final pose */
    Metric::ValueSequenceBase<float>* mDiffTranslation;
    /* Absolute difference between the initial angle and the final angle */
    Metric::ValueSequenceBase<float>* mDiffRotation;
    /* Total number of the iterations */
    Metric::ValueSequenceBase<int>*   mNumOfIterations;
    /* Total number of the step size updates */
    Metric::ValueSequenceBase<int>*   mNumOfRefinements;
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
    using ScanMatcherMetrics = ScanMatcherHillClimbingMetrics;

    /* Constructor */
    ScanMatcherHillClimbing(const std::string& scanMatcherName,
                            LikelihoodFuncPtr&& likelihoodFunc,
                            const double linearDelta,
                            const double angularDelta,
                            const int maxIterations,
                            const int numOfRefinements);

    /* Destructor */
    ~ScanMatcherHillClimbing() = default;

    /* Optimize particle poses by scan matching */
    ScanMatchingResultVector OptimizePose(
        const ScanMatchingQueryVector& queries,
        const Sensor::ScanDataPtr<double>& scanData) override;

private:
    /* Optimize particle pose by scan matching methods */
    ScanMatchingResult OptimizePoseCore(
        const std::size_t particleIdx,
        const ScanMatchingQuery& query,
        const Sensor::ScanDataPtr<double>& scanData);

private:
    /* Likelihood function to compute the observation likelihood */
    LikelihoodFuncPtr  mLikelihoodFunc;
    /* Initial step of the linear components (x and y) */
    const double       mLinearDelta;
    /* Initial step of the angular component (theta) */
    const double       mAngularDelta;
    /* Maximum number of the iterations */
    const int          mMaxIterations;
    /* Maximum number of the step parameter updates */
    const int          mNumOfRefinements;
    /* Metrics information */
    ScanMatcherMetrics mMetrics;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_MAPPING_SCAN_MATCHER_HILL_CLIMBING_HPP */
