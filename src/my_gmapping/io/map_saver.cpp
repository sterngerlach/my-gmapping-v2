
/* map_saver.cpp */

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "my_gmapping/io/map_saver.hpp"

/* Declare namespaces for convenience */
namespace gil = boost::gil;
namespace pt = boost::property_tree;

namespace MyGMapping {
namespace IO {

/* Save the specified map */
bool MapSaver::SaveMap(
    const Mapping::GridMapType& gridMap,
    const std::vector<Mapping::TimeStampedPose>& trajectory,
    const std::string& fileName,
    bool drawTrajectory) const
{
    /* Compute the actual map size */
    Point2D<int> patchIdxMin;
    Point2D<int> patchIdxMax;
    Point2D<int> gridCellIdxMin;
    Point2D<int> gridCellIdxMax;
    Point2D<int> mapSizeInPatches;
    Point2D<int> mapSizeInGridCells;
    gridMap.ComputeActualMapSize(patchIdxMin, patchIdxMax,
                                 gridCellIdxMin, gridCellIdxMax,
                                 mapSizeInPatches, mapSizeInGridCells);

    /* Initialize the grid map image */
    gil::rgb8_image_t mapImage {
        mapSizeInGridCells.mX, mapSizeInGridCells.mY };
    const gil::rgb8_view_t& mapImageView = gil::view(mapImage);
    gil::fill_pixels(mapImageView, gil::rgb8_pixel_t(192, 192, 192));

    /* Draw the grid cells in the map */
    this->DrawMap(gridMap, mapImageView, patchIdxMin, mapSizeInPatches);

    /* Draw the trajectory of the particle */
    if (drawTrajectory)
        this->DrawTrajectory(gridMap, trajectory, mapImageView,
                             gridCellIdxMin, mapSizeInGridCells);

    /* Save the map as PNG image
     * Image should be flipped upside down */
    try {
        const std::string pngFileName = fileName + ".png";
#if BOOST_VERSION <= 106700
        gil::png_write_view(pngFileName,
                            gil::flipped_up_down_view(mapImageView));
#else
        gil::write_view(pngFileName,
                        gil::flipped_up_down_view(mapImageView),
                        gil::png_tag());
#endif
    } catch (const std::ios_base::failure& e) {
        std::cerr << "std::ios_base::failure occurred: " << e.what() << ' '
                  << "(Error code: " << e.code() << ")" << std::endl;
        return false;
    }

    /* Save the map metadata as JSON format */
    try {
        const std::string metadataFileName = fileName + ".json";

        /* GridMap<T>::CellIndexToMapCoordinate(x, y) returns
         * the minimum position of the specified grid cell (x, y) */
        const Point2D<double> bottomLeft =
            gridMap.CellIndexToMapCoordinate(gridCellIdxMin);
        const Point2D<double> topRight =
            gridMap.CellIndexToMapCoordinate(gridCellIdxMax);

        /* Save the map metadata */
        this->SaveMapMetadata(gridMap.CellSize(), gridMap.PatchSize(),
                              mapSizeInPatches, mapSizeInGridCells,
                              bottomLeft, topRight, metadataFileName);
    } catch (const pt::json_parser_error& e) {
        std::cerr << "boost::property_tree::json_parser_error occurred: "
                  << e.what() << ' '
                  << "(" << e.filename() << ", Line " << e.line() << ")"
                  << std::endl;
        return false;
    }

    /* Save the robot trajectory */
    try {
        const std::string trajectoryFileName = fileName + ".trajectory";
        this->SaveTrajectory(trajectory, trajectoryFileName);
    } catch (const std::ios_base::failure& e) {
        std::cerr << "std::ios_base::failure occurred: " << e.what() << ' '
                  << "(Error code: " << e.code() << ")" << std::endl;
        return false;
    }

    return true;
}

/* Draw the grid cells to the image */
void MapSaver::DrawMap(const Mapping::GridMapType& gridMap,
                       const boost::gil::rgb8_view_t& mapImageView,
                       const Point2D<int>& patchIdxMin,
                       const Point2D<int>& mapSizeInPatches) const
{
    const double unknownVal = gridMap.UnknownValue();

    for (int y = 0; y < mapSizeInPatches.mY; ++y) {
        for (int x = 0; x < mapSizeInPatches.mX; ++x) {
            /* Retrieve the current map patch */
            const auto& patch = gridMap.PatchAt(
                patchIdxMin.mX + x, patchIdxMin.mY + y);

            if (!patch.IsAllocated())
                continue;

            /* Draw the grid cells in the patch */
            for (int yy = 0; yy < gridMap.PatchSize(); ++yy) {
                for (int xx = 0; xx < gridMap.PatchSize(); ++xx) {
                    const double gridCellVal = patch.At(xx, yy).Value();

                    /* If the occupancy probability value is less than or
                     * equal to zero, then the grid cell is not yet observed
                     * and is in unknown state (Cell::Unknown is zero) */
                    if (gridCellVal == unknownVal ||
                        gridCellVal < 0.0 || gridCellVal > 1.0)
                        continue;

                    const std::ptrdiff_t idxX = static_cast<std::ptrdiff_t>(
                        x * gridMap.PatchSize() + xx);
                    const std::ptrdiff_t idxY = static_cast<std::ptrdiff_t>(
                        y * gridMap.PatchSize() + yy);
                    const std::uint8_t grayScale = static_cast<std::uint8_t>(
                        (1.0 - gridCellVal) * 255.0);

                    mapImageView(idxX, idxY) =
                        gil::rgb8_pixel_t(grayScale, grayScale, grayScale);
                }
            }
        }
    }
}

/* Draw the trajectory of the particle to the image */
void MapSaver::DrawTrajectory(
    const Mapping::GridMapType& gridMap,
    const std::vector<Mapping::TimeStampedPose>& trajectory,
    const boost::gil::rgb8_view_t& mapImageView,
    const Point2D<int>& gridCellIdxMin,
    const Point2D<int>& mapSizeInGridCells) const
{
    /* Trajectory should contain at least one node */
    assert(!trajectory.empty());

    Point2D<int> prevGridCellIdx = gridMap.MapCoordinateToCellIndex(
        trajectory.front().mPose.mX, trajectory.front().mPose.mY);

    const Point2D<int> gridCellIdxMax {
        gridCellIdxMin.mX + mapSizeInGridCells.mX,
        gridCellIdxMin.mY + mapSizeInGridCells.mY };

    for (const auto& stampedPose : trajectory) {
        const Point2D<int> gridCellIdx = gridMap.MapCoordinateToCellIndex(
            stampedPose.mPose.mX, stampedPose.mPose.mY);
        std::vector<Point2D<int>> lineIndices;
        Bresenham(prevGridCellIdx, gridCellIdx, lineIndices);

        for (const auto& interpolatedIdx : lineIndices) {
            if (interpolatedIdx.mX < gridCellIdxMin.mX ||
                interpolatedIdx.mX >= gridCellIdxMax.mX - 1 ||
                interpolatedIdx.mY < gridCellIdxMin.mY ||
                interpolatedIdx.mY >= gridCellIdxMax.mY - 1)
                continue;

            const int x = interpolatedIdx.mX - gridCellIdxMin.mX;
            const int y = interpolatedIdx.mY - gridCellIdxMin.mY;
            gil::fill_pixels(gil::subimage_view(mapImageView, x, y, 2, 2),
                             gil::rgb8_pixel_t(255, 0, 0));
        }

        prevGridCellIdx = gridCellIdx;
    }
}

/* Save the robot trajectory */
void MapSaver::SaveTrajectory(
    const std::vector<Mapping::TimeStampedPose>& trajectory,
    const std::string& fileName) const
{
    /* Open the file (throw exception if failed) */
    std::ofstream outFile;
    outFile.exceptions(std::ios_base::failbit | std::ios_base::badbit);
    outFile.open(fileName);

    /* Set the floating-point number format */
    outFile << std::fixed << std::setprecision(6);

    /* Write the robot trajectory with timestamp */
    for (const auto& stampedPose : trajectory) {
        outFile << stampedPose.mTimeStamp << ' '
                << stampedPose.mPose.mX << ' '
                << stampedPose.mPose.mY << ' '
                << stampedPose.mPose.mTheta << '\n';
    }

    /* Close the file */
    outFile.close();

    return;
}

/* Save the map metadata as JSON format */
void MapSaver::SaveMapMetadata(
    const double mapResolution,
    const int patchSize,
    const Point2D<int>& mapSizeInPatches,
    const Point2D<int>& mapSizeInGridCells,
    const Point2D<double>& bottomLeft,
    const Point2D<double>& topRight,
    const std::string& fileName) const
{
    pt::ptree jsonTree;

    /* Write the map metadata */
    jsonTree.put("Map.Resolution", mapResolution);
    jsonTree.put("Map.PatchSize", patchSize);
    jsonTree.put("Map.WidthInPatches", mapSizeInPatches.mX);
    jsonTree.put("Map.HeightInPatches", mapSizeInPatches.mY);
    jsonTree.put("Map.WidthInGridCells", mapSizeInGridCells.mX);
    jsonTree.put("Map.HeightInGridCells", mapSizeInGridCells.mY);

    jsonTree.put("Map.BottomLeft.X", bottomLeft.mX);
    jsonTree.put("Map.BottomLeft.Y", bottomLeft.mY);
    jsonTree.put("Map.TopRight.X", topRight.mX);
    jsonTree.put("Map.TopRight.Y", topRight.mY);

    /* Save the json to the file */
    pt::write_json(fileName, jsonTree);

    return;
}

} /* namespace IO */
} /* namespace MyGMapping */
