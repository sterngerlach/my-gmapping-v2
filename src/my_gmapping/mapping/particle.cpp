
/* particle.cpp */

#include "my_gmapping/mapping/particle.hpp"

namespace MyGMapping {
namespace Mapping {

/* Constructor */
Particle::Particle(const GridMapType& gridMap,
                   const GridMapType& latestMap,
                   double weight,
                   const RobotPose2D<double>& pose,
                   const std::shared_ptr<TrajectoryNode>& node) :
    mMap(gridMap),
    mLatestMap(latestMap),
    mWeight(weight),
    mPose(pose),
    mNode(node)
{
}

/* Copy constructor */
Particle::Particle(const Particle& other) :
    mMap(other.mMap),
    mLatestMap(other.mLatestMap),
    mWeight(other.mWeight),
    mPose(other.mPose),
    mNode(other.mNode)
{
}

/* Copy assignment operator */
Particle& Particle::operator=(const Particle& other)
{
    if (this == &other)
        return *this;

    this->mMap = other.mMap;
    this->mLatestMap = other.mLatestMap;
    this->mWeight = other.mWeight;
    this->mPose = other.mPose;
    this->mNode = other.mNode;

    return *this;
}

/* Move constructor */
Particle::Particle(Particle&& other) :
    mMap(std::move(other.mMap)),
    mLatestMap(std::move(other.mLatestMap)),
    mWeight(other.mWeight),
    mPose(std::move(other.mPose)),
    mNode(std::move(other.mNode))
{
}

/* Move assignment operator */
Particle& Particle::operator=(Particle&& other)
{
    if (this == &other)
        return *this;

    this->mMap = std::move(other.mMap);
    this->mLatestMap = std::move(other.mLatestMap);
    this->mWeight = other.mWeight;
    this->mPose = std::move(other.mPose);
    this->mNode = std::move(other.mNode);

    return *this;
}

} /* namespace Mapping */
} /* namespace MyGMapping */
