
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

    /* Optimize particle pose by scan matching methods */
    virtual void OptimizePose(
        const std::size_t numOfParticles,
        const std::vector<const GridMap*>& particleMaps,
        const Sensor::ScanDataPtr<double>& scanData,
        const std::vector<RobotPose2D<double>>& initialPoses,
        std::vector<RobotPose2D<double>>& estimatedPoses,
        std::vector<double>& likelihoodValues) = 0;

protected:
    /* Name of this scan matcher */
    const std::string mName;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_SCAN_MATCHER_HPP */
