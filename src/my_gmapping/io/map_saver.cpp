
/* map_saver.cpp */

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "my_gmapping/io/map_saver.hpp"

#include "my_gmapping/bresenham.hpp"

/* Declare namespaces for convenience */
namespace gil = boost::gil;
namespace pt = boost::property_tree;

namespace MyGMapping {
namespace IO {

/* Save the specified map */
bool MapSaver::SaveMap(const Mapping::GridMap& gridMap,
                       const Trajectory& trajectory,
                       const std::string& fileName,
                       const bool drawTrajectory) const
{
    /* Compute the size of the cropped grid map */
    const BoundingBox<int> croppedBox = gridMap.CroppedBoundingBox();

    /* Initialize the grid map image */
    gil::rgb8_image_t mapImage { croppedBox.Width(), croppedBox.Height() };
    const gil::rgb8_view_t& mapImageView = gil::view(mapImage);
    gil::fill_pixels(mapImageView, gil::rgb8_pixel_t(192, 192, 192));

    /* Draw the grid cells in the map */
    this->DrawMap(mapImageView, gridMap, croppedBox);

    /* Draw the trajectory of the particle */
    if (drawTrajectory)
        this->DrawTrajectory(mapImageView, gridMap, trajectory, croppedBox);

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
        this->SaveMapMetadata(gridMap, croppedBox, metadataFileName);
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
void MapSaver::DrawMap(const boost::gil::rgb8_view_t& mapImageView,
                       const Mapping::GridMap& gridMap,
                       const BoundingBox<int>& boundingBox) const
{
    /* Check that the grid map has the same size as the image */
    Assert(mapImageView.width() == boundingBox.Width());
    Assert(mapImageView.height() == boundingBox.Height());

    const auto toGrayScale = [](const std::uint8_t value) {
        return gil::rgb8_pixel_t(value, value, value); };
    const auto unknownProb = gridMap.UnknownProbability();

    /* Draw the grid cells to the image */
    for (int row = boundingBox.mMin.mY; row < boundingBox.mMax.mY; ++row) {
        for (int col = boundingBox.mMin.mX; col < boundingBox.mMax.mX; ++col) {
            /* Get the probability value */
            const auto prob = gridMap.ProbabilityOr(row, col, unknownProb);

            /* Skip if the grid cell is unknown */
            if (prob == unknownProb)
                continue;

            /* Convert the probability to the pixel intensity */
            const auto imageRow = static_cast<std::ptrdiff_t>(
                row - boundingBox.mMin.mY);
            const auto imageCol = static_cast<std::ptrdiff_t>(
                col - boundingBox.mMin.mX);
            const auto pixelValue = static_cast<std::uint8_t>(
                (1.0 - prob) * 255.0);

            /* Set the pixel intensity */
            mapImageView(imageCol, imageRow) = toGrayScale(pixelValue);
        }
    }
}

/* Draw the trajectory of the particle to the image */
void MapSaver::DrawTrajectory(const boost::gil::rgb8_view_t& mapImageView,
                              const Mapping::GridMap& gridMap,
                              const Trajectory& trajectory,
                              const BoundingBox<int>& boundingBox) const
{
    /* Trajectory should contain at least one node */
    if (trajectory.size() < 2)
        return;

    Point2D<int> prevIdx = gridMap.PositionToIndex(
        trajectory.front().mPose.mX, trajectory.front().mPose.mY);

    for (const auto& stampedPose : trajectory) {
        const Point2D<int> currentIdx = gridMap.PositionToIndex(
            stampedPose.mPose.mX, stampedPose.mPose.mY);
        std::vector<Point2D<int>> lineIndices;
        Bresenham(prevIdx, currentIdx, lineIndices);

        for (const auto& interpolatedIdx : lineIndices) {
            if (interpolatedIdx.mX < boundingBox.mMin.mX ||
                interpolatedIdx.mX > boundingBox.mMax.mX - 1 ||
                interpolatedIdx.mY < boundingBox.mMin.mY ||
                interpolatedIdx.mY > boundingBox.mMax.mY - 1)
                continue;

            const int x = interpolatedIdx.mX - boundingBox.mMin.mX;
            const int y = interpolatedIdx.mY - boundingBox.mMin.mY;
            gil::fill_pixels(gil::subimage_view(mapImageView, x, y, 2, 2),
                             gil::rgb8_pixel_t(255, 0, 0));
        }

        prevIdx = currentIdx;
    }
}

/* Save the robot trajectory */
void MapSaver::SaveTrajectory(const Trajectory& trajectory,
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
void MapSaver::SaveMapMetadata(const Mapping::GridMap& gridMap,
                               const BoundingBox<int>& boundingBox,
                               const std::string& fileName) const
{
    pt::ptree jsonMetadata;

    /* Write the map metadata */
    jsonMetadata.put("Resolution", gridMap.Resolution());
    jsonMetadata.put("Log2BlockSize", gridMap.Log2BlockSize());
    jsonMetadata.put("BlockSize", gridMap.BlockSize());

    const int rows = boundingBox.Height();
    const int cols = boundingBox.Width();
    const int blockRows = (rows + gridMap.BlockSize() - 1)
                          >> gridMap.Log2BlockSize();
    const int blockCols = (cols + gridMap.BlockSize() - 1)
                          >> gridMap.Log2BlockSize();
    const double height = gridMap.Resolution() * rows;
    const double width = gridMap.Resolution() * cols;

    jsonMetadata.put("Rows", rows);
    jsonMetadata.put("Cols", cols);
    jsonMetadata.put("BlockRows", blockRows);
    jsonMetadata.put("BlockCols", blockCols);
    jsonMetadata.put("Height", height);
    jsonMetadata.put("Width", width);

    const auto posMin = gridMap.IndexToPosition(
        boundingBox.mMin.mY, boundingBox.mMin.mX);
    const auto posMax = gridMap.IndexToPosition(
        boundingBox.mMax.mY, boundingBox.mMax.mX);

    jsonMetadata.put("PositionMin.X", posMin.mX);
    jsonMetadata.put("PositionMin.Y", posMin.mY);
    jsonMetadata.put("PositionMax.X", posMax.mX);
    jsonMetadata.put("PositionMax.Y", posMax.mY);

    /* Save the json to the file */
    pt::write_json(fileName, jsonMetadata);

    return;
}

} /* namespace IO */
} /* namespace MyGMapping */
