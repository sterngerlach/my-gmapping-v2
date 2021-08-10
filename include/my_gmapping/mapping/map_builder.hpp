
/* map_builder.hpp */

#ifndef MY_GMAPPING_MAP_BUILDER_HPP
#define MY_GMAPPING_MAP_BUILDER_HPP

#include <deque>
#include <vector>

#include "my_gmapping/grid_map/grid_map.hpp"
#include "my_gmapping/grid_map/grid_map_binary_bayes_grid_cell.hpp"
#include "my_gmapping/grid_map/grid_map_cell.hpp"
#include "my_gmapping/point.hpp"
#include "my_gmapping/pose.hpp"
#include "my_gmapping/util.hpp"
#include "my_gmapping/mapping/grid_map_types.hpp"
#include "my_gmapping/sensor/sensor_data.hpp"

namespace MyGMapping {
namespace Mapping {

class MapBuilder
{
public:
    /* Type definitions for convenience */
    using RobotPoseDeque = std::deque<RobotPose2D<double>>;
    using ScanDataDeque = std::deque<Sensor::ScanDataPtr<double>>;

    /* Constructor */
    MapBuilder(const double maxUsableRange,
               const double minUsableRange,
               const double probHit,
               const double probMiss);
    /* Destructor */
    ~MapBuilder() = default;

    /* Copy constructor (disabled) */
    MapBuilder(const MapBuilder&) = delete;
    /* Copy assignment operator (disabled) */
    MapBuilder& operator=(const MapBuilder&) = delete;
    /* Move constructor */
    MapBuilder(MapBuilder&&) = default;
    /* Move assignment operator */
    MapBuilder& operator=(MapBuilder&&) = default;

    /* Integrate the scan data to the grid map */
    /* Update the particle map using the latest scan data */
    void UpdateGridMap(GridMap& gridMap,
                       const RobotPose2D<double>& currentPose,
                       const Sensor::ScanDataPtr<double>& scanData) const;

    /* Update the particle map with the multiple latest scans */
    void UpdateLatestMap(GridMap& latestMap,
                         const RobotPoseDeque& latestPoses,
                         const ScanDataDeque& latestScanData);

    /* Compute the indices of the missed cells using Bresenham algorithm */
    void ComputeMissedCellIndices(
        const Point2D<int>& startCellIdx,
        const Point2D<int>& endCellIdx,
        std::vector<Point2D<int>>& gridCellIndices) const;

    /* Compute the indices of the missed cells using the Bresenham algorithm
     * at the subpixel accuracy */
    void ComputeMissedIndicesScaled(
        const Point2D<int>& scaledStartIdx,
        const Point2D<int>& scaledEndIdx,
        const int subpixelScale,
        std::vector<Point2D<int>>& missedIndices) const;

    /* Get the subpixel scale used when computing missed grid cell indices */
    inline int RayCastingSubpixelScale() const { return SubpixelScale; }
    /* Get the minimum scan range considered valid when creating a grid map */
    inline double MinRange() const { return this->mMinUsableRange; }
    /* Get the maximum scan range considered valid when creating a grid map */
    inline double MaxRange() const { return this->mMaxUsableRange; }
    /* Get the probability for hit grid cells */
    inline double ProbabilityHit() const { return this->mProbHit; }
    /* Get the probability for missed grid cells */
    inline double ProbabilityMiss() const { return this->mProbMiss; }
    /* Get the odds for hit grid cells */
    inline double OddsHit() const { return this->mOddsHit; }
    /* Get the odds for missed grid cells */
    inline double OddsMiss() const { return this->mOddsMiss; }

private:
    /* Subpixel scale for computing the missed grid cell indices */
    static constexpr int SubpixelScale = 100;

private:
    /* Maximum scan range considered to be valid */
    const double mMaxUsableRange;
    /* Minimum scan range considered to be valid */
    const double mMinUsableRange;
    /* Probability for hit grid cells */
    const double mProbHit;
    /* Probability for missed grid cells */
    const double mProbMiss;
    /* Odds for hit grid cells */
    const double mOddsHit;
    /* Odds for missed grid cells */
    const double mOddsMiss;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_MAP_BUILDER_HPP */
