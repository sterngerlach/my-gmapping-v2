
/* map_saver.hpp */

#ifndef MY_GMAPPING_MAP_SAVER_HPP
#define MY_GMAPPING_MAP_SAVER_HPP

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

/* Definitions to prevent compile errors
 * int_p_NULL is removed in libpng 1.4 */
#define png_infopp_NULL     (png_infopp)NULL
#define int_p_NULL          (int*)NULL
#define png_bytep_NULL      (png_bytep)NULL

#include <boost/version.hpp>

#if BOOST_VERSION <= 106700
#include <boost/gil/gil_all.hpp>
#include <boost/gil/extension/io/png_io.hpp>
#elif BOOST_VERSION <= 106800
#include <boost/gil/gil_all.hpp>
#include <boost/gil/extension/io/png.hpp>
#else
#include <boost/gil.hpp>
#include <boost/gil/extension/io/png.hpp>
#endif

#include "my_gmapping/point.hpp"
#include "my_gmapping/pose.hpp"
#include "my_gmapping/util.hpp"
#include "my_gmapping/grid_map/grid_map.hpp"
#include "my_gmapping/grid_map/grid_map_patch.hpp"
#include "my_gmapping/mapping/grid_map_builder.hpp"
#include "my_gmapping/mapping/trajectory_node.hpp"

namespace MyGMapping {
namespace IO {

class MapSaver
{
public:
    /* Constructor */
    MapSaver() = default;

    /* Destructor */
    ~MapSaver() = default;

    /* Save the specified map */
    bool SaveMap(const Mapping::GridMapType& gridMap,
                 const std::vector<Mapping::TimeStampedPose>& trajectory,
                 const std::string& fileName,
                 bool drawTrajectory) const;

private:
    /* Draw the grid cells to the image */
    void DrawMap(const Mapping::GridMapType& gridMap,
                 const boost::gil::rgb8_view_t& mapImageView,
                 const Point2D<int>& patchIdxMin,
                 const Point2D<int>& mapSizeInPatches) const;

    /* Draw the trajectory of the particle to the image */
    void DrawTrajectory(const Mapping::GridMapType& gridMap,
                        const std::vector<Mapping::TimeStampedPose>& trajectory,
                        const boost::gil::rgb8_view_t& mapImageView,
                        const Point2D<int>& gridCellIdxMin,
                        const Point2D<int>& mapSizeInGridCells) const;

    /* Save the robot trajectory */
    void SaveTrajectory(
        const std::vector<Mapping::TimeStampedPose>& trajectory,
        const std::string& fileName) const;

    /* Save the map metadata as JSON format */
    void SaveMapMetadata(const double mapResolution,
                         const int patchSize,
                         const Point2D<int>& mapSizeInPatches,
                         const Point2D<int>& mapSizeInGridCells,
                         const Point2D<double>& bottomLeft,
                         const Point2D<double>& topRight,
                         const std::string& fileName) const;
};

} /* namespace IO */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_MAP_SAVER_HPP */
