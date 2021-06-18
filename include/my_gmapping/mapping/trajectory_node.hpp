
/* trajectory_node.hpp */

#ifndef MY_GMAPPING_TRAJECTORY_NODE_HPP
#define MY_GMAPPING_TRAJECTORY_NODE_HPP

#include <memory>

#include "my_gmapping/pose.hpp"

namespace MyGMapping {
namespace Mapping {

/*
 * Robot pose with timestamp
 */
struct TimeStampedPose
{
    /* Constructor */
    TimeStampedPose(const RobotPose2D<double>& pose,
                    const double timeStamp) :
        mPose(pose), mTimeStamp(timeStamp) { }

    RobotPose2D<double> mPose;
    double              mTimeStamp;
};

/*
 * TrajectoryNode holds the particle trajectory as the tree structure
 */
class TrajectoryNode final
{
public:
    /* Constructor */
    TrajectoryNode(const std::shared_ptr<TrajectoryNode>& pParent,
                   const RobotPose2D<double>& pose,
                   const double timeStamp) :
        mParent(pParent), mStampedPose(pose, timeStamp) { }
    /* Destructor */
    ~TrajectoryNode() = default;

    /* Get the parent of the trajectory node */
    inline const std::shared_ptr<TrajectoryNode>& Parent() const
    { return this->mParent; }

    /* Get the pose with time stamp */
    inline const TimeStampedPose& StampedPose() const
    { return this->mStampedPose; }

    /* Get the pose of the trajectory node */
    inline const RobotPose2D<double>& Pose() const
    { return this->mStampedPose.mPose; }

    /* Get the time stamp */
    inline double TimeStamp() const
    { return this->mStampedPose.mTimeStamp; }

private:
    std::shared_ptr<TrajectoryNode> mParent;
    const TimeStampedPose           mStampedPose;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_TRAJECTORY_NODE_HPP */
