
/* particle.hpp */

#ifndef MY_GMAPPING_PARTICLE_HPP
#define MY_GMAPPING_PARTICLE_HPP

#include "my_gmapping/grid_map/grid_map.hpp"
#include "my_gmapping/grid_map/grid_map_cell.hpp"
#include "my_gmapping/grid_map/grid_map_binary_bayes_grid_cell.hpp"
#include "my_gmapping/pose.hpp"
#include "my_gmapping/mapping/grid_map_types.hpp"
#include "my_gmapping/mapping/trajectory_node.hpp"

namespace MyGMapping {
namespace Mapping {

class Particle final
{
public:
    /* Type definitions */
    using TrajectoryNodePtr = std::shared_ptr<TrajectoryNode>;

    /* Constructor */
    Particle(const GridMap& gridMap,
             const GridMap& latestMap,
             const double weight,
             const RobotPose2D<double>& pose,
             const TrajectoryNodePtr& node) :
        mMap(gridMap), mLatestMap(latestMap),
        mWeight(weight), mPose(pose), mNode(node) { }
    /* Destructor */
    ~Particle() = default;

    /* Copy constructor */
    Particle(const Particle& other) = default;
    /* Copy assignment operator */
    Particle& operator=(const Particle& other) = default;
    /* Move constructor */
    Particle(Particle&& other) = default;
    /* Move assignment operator */
    Particle& operator=(Particle&& other) = default;

    /* Get the occupancy grid map of this particle */
    inline GridMap& Map() { return this->mMap; }
    /* Get the occupancy grid map of this particle */
    inline const GridMap& Map() const { return this->mMap; }

    /* Get the latest grid map of this particle */
    inline GridMap& LatestMap() { return this->mLatestMap; }
    /* Get the latest grid map of this particle */
    inline const GridMap& LatestMap() const { return this->mLatestMap; }

    /* Get the particle weight */
    inline double Weight() const { return this->mWeight; }
    /* Set the particle weight */
    inline void SetWeight(const double weight) { this->mWeight = weight; }

    /* Get the trajectory node */
    inline TrajectoryNodePtr& Node() { return this->mNode; }
    /* Get the trajectory node */
    inline const TrajectoryNodePtr& Node() const { return this->mNode; }

    /* Get the current particle pose */
    inline RobotPose2D<double>& Pose() { return this->mPose; }
    /* Get the current particle pose */
    inline const RobotPose2D<double>& Pose() const { return this->mPose; }

private:
    /* Occupancy grid map associated to the particle */
    GridMap             mMap;
    /* Occupancy grid map constructed from the latest scans */
    GridMap             mLatestMap;
    /* Weight of the particle */
    double              mWeight;
    /* Current particle pose */
    RobotPose2D<double> mPose;
    /* Pointer to the trajectory node */
    TrajectoryNodePtr   mNode;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_PARTICLE_HPP */
