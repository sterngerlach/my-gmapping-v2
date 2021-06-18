
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

private:
    /* Subpixel scale for computing the missed grid cell indices */
    static constexpr int SubpixelScale = 100;

private:
    /* Maximum scan range considered to be valid */
    const double mMaxUsableRange;
    /* Minimum scan range considered to be valid */
    const double mMinUsableRange;
    /* Occupancy probability value for hit grid cell
     * Used for updating the value with Binary Bayes Filter */
    const double mProbHit;
    /* Occupancy probability value for missed grid cell
     * Used for updating the value with Binary Bayes Filter */
    const double mProbMiss;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_MAP_BUILDER_HPP */
