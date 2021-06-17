
/* scan_interpolator.hpp */

#ifndef MY_GMAPPING_MAPPING_SCAN_INTERPOLATOR_HPP
#define MY_GMAPPING_MAPPING_SCAN_INTERPOLATOR_HPP

#include <cassert>
#include <memory>
#include <vector>

#include "my_gmapping/point.hpp"
#include "my_gmapping/pose.hpp"
#include "my_gmapping/util.hpp"
#include "my_gmapping/sensor/sensor_data.hpp"

namespace MyGMapping {
namespace Mapping {

/* Type definition for convenience */
class ScanInterpolator;
using ScanInterpolatorPtr = std::unique_ptr<ScanInterpolator>;

class ScanInterpolator final
{
public:
    /* Constructor */
    ScanInterpolator(const double distScans,
                     const double distThresholdEmpty) :
        mDistScans(distScans),
        mDistThresholdEmpty(distThresholdEmpty) { }

    /* Destructor */
    ~ScanInterpolator() = default;

    /* Interpolate scan data */
    Sensor::ScanDataPtr<double> Interpolate(
        const Sensor::ScanDataPtr<double>& scanData) const;

private:
    /* Distance between two scan points */
    const double mDistScans;
    /* Distance threshold for empty space */
    const double mDistThresholdEmpty;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_MAPPING_SCAN_INTERPOLATOR_HPP */
