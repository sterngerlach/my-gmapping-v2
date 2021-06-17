
/* particle.hpp */

#ifndef MY_GMAPPING_PARTICLE_HPP
#define MY_GMAPPING_PARTICLE_HPP

#include "my_gmapping/pose.hpp"
#include "my_gmapping/mapping/grid_map_types.hpp"
#include "my_gmapping/mapping/trajectory_node.hpp"

namespace MyGMapping {
namespace Mapping {

struct Particle final
{
    /* Constructor */
    Particle(const GridMap& gridMap,
             const GridMap& latestMap,
             const double weight,
             const RobotPose2D<double>& pose,
             const std::shared_ptr<TrajectoryNode>& node) :
        mMap(gridMap), mLatestMap(latestMap),
        mWeight(weight), mPose(pose), mNode(node) { }
    /* Destructor */
    ~Particle() = default;

    /* Occupancy grid map associated to the particle */
    GridMap mMap;
    /* Occupancy grid map constructed from the latest scans */
    GridMap mLatestMap;
    /* Weight of the particle */
    double  mWeight;
    /* Current particle pose */
    RobotPose2D<double> mPose;
    /* Pointer to the trajectory node */
    std::shared_ptr<TrajectoryNode> mNode;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_PARTICLE_HPP */
