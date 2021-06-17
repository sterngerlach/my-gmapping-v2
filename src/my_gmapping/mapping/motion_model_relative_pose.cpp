
/* motion_model_relative_pose.cpp */

#include "my_gmapping/mapping/motion_model_relative_pose.hpp"

#include "my_gmapping/util.hpp"

namespace MyGMapping {
namespace Mapping {

/* Sample value from normal distribution */
double MotionModelRelativePose::StandardNormalDist() const
{
    /* Standard normal distribution used for motion sampling */
    static std::normal_distribution stdNormalDist { 0.0, 1.0 };
    return stdNormalDist(this->mRandEngine);
}

/* Sample new pose using the previous pose and the latest control */
RobotPose2D<double> MotionModelRelativePose::SamplePose(
    const RobotPose2D<double>& prevPose,
    const RobotPose2D<double>& motionControl)
{
    const double sigmaXY = this->mSigmaLinear * 0.25;

    const double sigmaX = std::sqrt(
        this->mSigmaLinear * std::fabs(motionControl.mX) +
        this->mSigmaAngularToLinear * std::fabs(motionControl.mTheta) +
        sigmaXY * std::fabs(motionControl.mY));
    const double sigmaY = std::sqrt(
        this->mSigmaLinear * std::fabs(motionControl.mY) +
        this->mSigmaAngularToLinear * std::fabs(motionControl.mTheta) +
        sigmaXY * std::fabs(motionControl.mX));
    const double sigmaTheta = std::sqrt(
        this->mSigmaAngular * std::fabs(motionControl.mTheta) +
        this->mSigmaLinearToAngular * std::sqrt(
            motionControl.mX * motionControl.mX +
            motionControl.mY * motionControl.mY));

    const RobotPose2D<double> noisedControl {
        motionControl.mX + this->StandardNormalDist() * sigmaX,
        motionControl.mY + this->StandardNormalDist() * sigmaY,
        motionControl.mTheta + this->StandardNormalDist() * sigmaTheta };

    return Compound(prevPose, noisedControl);
}

} /* namespace Mapping */
} /* namespace MyGMapping */
