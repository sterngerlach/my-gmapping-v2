
/* motion_model.hpp */

#ifndef MY_GMAPPING_MOTION_MODEL_HPP
#define MY_GMAPPING_MOTION_MODEL_HPP

#include <cmath>
#include <memory>

#include "my_gmapping/pose.hpp"

namespace MyGMapping {
namespace Mapping {

/* Type definitions for convenience */
class MotionModel;
using MotionModelPtr = std::unique_ptr<MotionModel>;

class MotionModel
{
public:
    /* Constructor */
    MotionModel() = default;
    /* Destructor */
    virtual ~MotionModel() = default;

    /* Copy constructor (disabled) */
    MotionModel(const MotionModel&) = delete;
    /* Copy assignment operator (disabled) */
    MotionModel& operator=(const MotionModel&) = delete;
    /* Move constructor (disabled) */
    MotionModel(MotionModel&&) = delete;
    /* Move assignment operator (disabled) */
    MotionModel& operator=(MotionModel&&) = delete;

    /* Sample new pose */
    virtual RobotPose2D<double> SamplePose(
        const RobotPose2D<double>& prevPose,
        const RobotPose2D<double>& motionControl) = 0;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_MOTION_MODEL_HPP */
