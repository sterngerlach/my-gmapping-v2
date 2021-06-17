
/* scan_matcher.hpp */

#ifndef MY_GMAPPING_SCAN_MATCHER_HPP
#define MY_GMAPPING_SCAN_MATCHER_HPP

#include <cassert>
#include <cmath>
#include <functional>
#include <memory>
#include <vector>

#include "my_gmapping/point.hpp"
#include "my_gmapping/pose.hpp"
#include "my_gmapping/mapping/grid_map_types.hpp"
#include "my_gmapping/sensor/sensor_data.hpp"

namespace MyGMapping {
namespace Mapping {

/* Type definition for convenience */
class ScanMatcher;
using ScanMatcherPtr = std::unique_ptr<ScanMatcher>;

/*
 * ScanMatchingQuery struct holds the necessary information for performing
 * scan-to-map scan matching, such as an initial robot pose in a world
 * coordinate frame and a grid map
 */
struct ScanMatchingQuery final
{
    /* Constructor */
    ScanMatchingQuery(const GridMap& pGridMap,
                      const RobotPose2D<double>& initialPose) :
        mGridMap(pGridMap),
        mInitialPose(initialPose) { }

    /* Destructor */
    ~ScanMatchingQuery() = default;

    /* Grid map */
    const GridMap&                    mGridMap;
    /* Initial robot pose in a world coordinate frame */
    const RobotPose2D<double>         mInitialPose;
};

/* Vector of the scan matching query */
using ScanMatchingQueryVector = std::vector<ScanMatchingQuery>;

/*
 * ScanMatchingResult struct holds the details of the scan matching result
 */
struct ScanMatchingResult
{
    /* Default constructor */
    ScanMatchingResult() = default;

    /* Constructor */
    ScanMatchingResult(const RobotPose2D<double>& initialPose,
                        const RobotPose2D<double>& estimatedPose,
                        const double normalizedLikelihood,
                        const double likelihood,
                        const double normalizedScore,
                        const double score) :
        mInitialPose(initialPose),
        mEstimatedPose(estimatedPose),
        mNormalizedLikelihood(normalizedLikelihood),
        mLikelihood(likelihood),
        mNormalizedScore(normalizedScore),
        mScore(score) { }

    /* Destructor */
    ~ScanMatchingResult() = default;

    /* Initial robot pose in a world coordinate frame */
    RobotPose2D<double> mInitialPose;
    /* Estimated robot pose in a world coordinate frame */
    RobotPose2D<double> mEstimatedPose;
    /* Normalized likelihood value */
    double              mNormalizedLikelihood;
    /* Likelihood value */
    double              mLikelihood;
    /* Normalized score value */
    double              mNormalizedScore;
    /* Score value */
    double              mScore;
};

/* Vector of the scan matching summary */
using ScanMatchingResultVector = std::vector<ScanMatchingResult>;

class ScanMatcher
{
public:
    /* Constructor */
    ScanMatcher(const std::string& scanMatcherName) :
        mName(scanMatcherName) { }
    /* Destructor */
    virtual ~ScanMatcher() = default;

    /* Retrieve the name of this scan matcher */
    inline const std::string& Name() const { return this->mName; }

    /* Optimize particle poses by scan matching */
    virtual ScanMatchingResultVector OptimizePose(
        const ScanMatchingQueryVector& queries,
        const Sensor::ScanDataPtr<double>& scanData) = 0;

protected:
    /* Name of this scan matcher */
    const std::string mName;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_SCAN_MATCHER_HPP */
