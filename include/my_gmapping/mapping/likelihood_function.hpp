
/* likelihood_function.hpp */

#ifndef MY_GMAPPING_LIKELIHOOD_FUNCTION_HPP
#define MY_GMAPPING_LIKELIHOOD_FUNCTION_HPP

#include <cmath>
#include <cstdlib>
#include <memory>

#include "my_gmapping/point.hpp"
#include "my_gmapping/pose.hpp"
#include "my_gmapping/mapping/grid_map_types.hpp"
#include "my_gmapping/sensor/sensor_data.hpp"

namespace MyGMapping {
namespace Mapping {

/* Type definition for convenience */
class LikelihoodFunction;
using LikelihoodFunctionPtr = std::unique_ptr<LikelihoodFunction>;

class LikelihoodFunction
{
public:
    /* Constructor */
    LikelihoodFunction() = default;
    /* Destructor */
    virtual ~LikelihoodFunction() = default;

    /* Copy constructor (disabled) */
    LikelihoodFunction(const LikelihoodFunction&) = delete;
    /* Copy assignment operator (disabled) */
    LikelihoodFunction& operator=(const LikelihoodFunction&) = delete;
    /* Move constructor (disabled) */
    LikelihoodFunction(LikelihoodFunction&&) = delete;
    /* Move assignment operator (disabled) */
    LikelihoodFunction& operator=(LikelihoodFunction&&) = delete;

    /* Calculate observation likelihood */
    virtual double Likelihood(
        const GridMapInterface& gridMap,
        const Sensor::ScanDataPtr<double>& scanData,
        const RobotPose2D<double>& sensorPose) = 0;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_LIKELIHOOD_FUNCTION_HPP */
