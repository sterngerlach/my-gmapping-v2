
/* motion_model_relative_pose.hpp */

#ifndef MY_GMAPPING_MAPPING_MOTION_MODEL_RELATIVE_POSE_HPP
#define MY_GMAPPING_MAPPING_MOTION_MODEL_RELATIVE_POSE_HPP

#include "my_gmapping/mapping/motion_model.hpp"

#include <random>

namespace MyGMapping {
namespace Mapping {

class MotionModelRelativePose final : public MotionModel
{
public:
    /* Constructor */
    MotionModelRelativePose(
        const double sigmaLinear,
        const double sigmaAngular,
        const double sigmaLinearToAngular,
        const double sigmaAngularToLinear) :
        MotionModel(),
        mSigmaLinear(sigmaLinear),
        mSigmaAngular(sigmaAngular),
        mSigmaLinearToAngular(sigmaLinearToAngular),
        mSigmaAngularToLinear(sigmaAngularToLinear),
        mRandEngine(std::random_device()()) { }

    /* Destructor */
    ~MotionModelRelativePose() = default;

    /* Sample new pose */
    RobotPose2D<double> SamplePose(
        const RobotPose2D<double>& prevPose,
        const RobotPose2D<double>& motionControl) override;

private:
    /* Sample value from standard normal distribution */
    double StandardNormalDist() const;

private:
    const double         mSigmaLinear;
    const double         mSigmaAngular;
    const double         mSigmaLinearToAngular;
    const double         mSigmaAngularToLinear;
    mutable std::mt19937 mRandEngine;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_MAPPING_MOTION_MODEL_RELATIVE_POSE_HPP */
