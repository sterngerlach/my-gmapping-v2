
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
    /* Constructor */
    Particle(const GridMapType& gridMap,
             const GridMapType& latestMap,
             double weight,
             const RobotPose2D<double>& pose,
             const std::shared_ptr<TrajectoryNode>& node);
    /* Destructor */
    ~Particle() = default;

    /* Copy constructor */
    Particle(const Particle& other);
    /* Copy assignment operator */
    Particle& operator=(const Particle& other);

    /* Move constructor */
    Particle(Particle&& other);
    /* Move assignment operator */
    Particle& operator=(Particle&& other);

    /* Get the occupancy grid map */
    inline GridMapType& Map() { return this->mMap; }
    inline const GridMapType& Map() const { return this->mMap; }

    /* Get the latest grid map */
    inline GridMapType& LatestMap() { return this->mLatestMap; }
    inline const GridMapType& LatestMap() const { return this->mLatestMap; }

    /* Get the particle weight */
    inline double Weight() const { return this->mWeight; }
    inline void SetWeight(double weight) { this->mWeight = weight; }

    /* Get the trajectory node */
    inline std::shared_ptr<TrajectoryNode>& Node()
    { return this->mNode; }
    inline const std::shared_ptr<TrajectoryNode>& Node() const
    { return this->mNode; }

    /* Get the current particle pose */
    inline RobotPose2D<double>& Pose() { return this->mPose; }
    inline const RobotPose2D<double>& Pose() const { return this->mPose; }

private:
    /* Occupancy grid map associated to the particle */
    GridMapType                     mMap;
    /* Occupancy grid map constructed from the latest scans */
    GridMapType                     mLatestMap;
    /* Weight of the particle */
    double                          mWeight;
    /* Current particle pose */
    RobotPose2D<double>             mPose;
    /* Pointer to the trajectory node */
    std::shared_ptr<TrajectoryNode> mNode;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_PARTICLE_HPP */
