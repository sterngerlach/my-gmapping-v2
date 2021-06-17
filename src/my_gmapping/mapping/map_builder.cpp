
/* map_builder.cpp */

#include "my_gmapping/mapping/map_builder.hpp"

namespace MyGMapping {
namespace Mapping {

/* Constructor */
MapBuilder::MapBuilder(const double maxUsableRange,
                       const double minUsableRange,
                       const double probHit,
                       const double probMiss) :
    mMaxUsableRange(maxUsableRange),
    mMinUsableRange(minUsableRange),
    mProbHit(probHit),
    mProbMiss(probMiss)
{
    XAssert(minUsableRange >= 0.0,
            "The minimum range of the laser scan that is considered valid "
            "must be larger than or equal to 0");
    XAssert(minUsableRange < maxUsableRange,
            "The maximum range of the laser scan that is considered valid "
            "must be larger than the minimum range");
    XAssert(probHit <= 1.0,
            "The probability value used for updating the hit grid cell "
            "with the Binary Bayes Filter must be less than or equal to 1.0");
    XAssert(probMiss >= 0.0,
            "The probability value used for updating the missed grid cell "
            "with the Binary Bayes Filter must be greater than or equal to 0.0");
    XAssert(probHit > probMiss,
            "The occupancy probability value for the hit grid cell must be "
            "larger than the one for the missed grid cell");
}

/* Integrate the scan data to the grid map */
/* Update the particle map using the latest scan data */
void MapBuilder::UpdateGridMap(
    GridMapType& gridMap,
    const RobotPose2D<double>& currentPose,
    const Sensor::ScanDataPtr<double>& scanData) const
{
    /* Vector for storing missed grid cell indices
     * Specified as static variable to reduce the performance loss
     * caused by memory allocations and deallocations */
    static std::vector<Point2D<int>> missedCellIndices;

    /* Vector for storing hit points in the world coordinate
     * Specified as static variable to reduce the performance loss
     * caused by memory allocations and deallocations */
    static std::vector<Point2D<double>> hitPoints;

    /* Compute the sensor pose from the particle pose
     * Do not calculate the cell index corresponding to the sensor pose here
     * since the map is later updated and anchor position is changed */
    const RobotPose2D<double> sensorPose =
        Compound(currentPose, scanData->RelativeSensorPose());

    /* Bottom-left corner position of the map */
    Point2D<double> bottomLeft { sensorPose.mX, sensorPose.mY };
    /* Top-right corner position of the map */
    Point2D<double> topRight { sensorPose.mX, sensorPose.mY };

    /* Store hit points in the world coordinate to reuse */
    hitPoints.clear();
    hitPoints.reserve(scanData->NumOfScans());

    /* Minimum range and maximum range */
    const double minRange = std::max(
        this->mMinUsableRange, scanData->MinRange());
    const double maxRange = std::min(
        this->mMaxUsableRange, scanData->MaxRange());

    /* Expand the map so that the particle map
     * contains all the scan points */
    for (std::size_t i = 0; i < scanData->NumOfScans(); ++i) {
        const double range = scanData->RangeAt(i);

        if (range >= maxRange || range <= minRange)
            continue;

        /* Compute the scan point */
        Point2D<double> hitPoint = scanData->HitPoint(sensorPose, i);

        /* Update the corner positions */
        bottomLeft.mX = std::min(bottomLeft.mX, hitPoint.mX);
        bottomLeft.mY = std::min(bottomLeft.mY, hitPoint.mY);
        topRight.mX = std::max(topRight.mX, hitPoint.mX);
        topRight.mY = std::max(topRight.mY, hitPoint.mY);

        /* Append the hit point */
        hitPoints.emplace_back(std::move(hitPoint));
    }

    /* Expand the map if necessary */
    gridMap.Expand(bottomLeft.mX, topRight.mX,
                   bottomLeft.mY, topRight.mY);

    /* Calculate the cell index corresponding to the sensor pose
     * since the map is updated and index might be changed */
    const Point2D<int> sensorCellIdx =
        gridMap.MapCoordinateToCellIndex(sensorPose.mX, sensorPose.mY);

    /* Integrate the scan into the particle map */
    const std::size_t numOfFilteredScans = hitPoints.size();

    for (std::size_t i = 0; i < numOfFilteredScans; ++i) {
        /* Compute the index of the hit cell */
        const Point2D<int> hitCellIdx =
            gridMap.MapCoordinateToCellIndex(hitPoints[i]);

        /* Compute the indices of the missed cells */
        this->ComputeMissedCellIndices(
            sensorCellIdx, hitCellIdx, missedCellIndices);

        /* Update the cell value */
        const std::size_t numOfMissedCells = missedCellIndices.size();

        for (std::size_t j = 0; j < numOfMissedCells; ++j)
            gridMap.Update(missedCellIndices[j], this->mProbMiss);

        gridMap.Update(hitCellIdx, this->mProbHit);
    }
}

/* Update the particle map with the multiple latest scans */
void MapBuilder::UpdateLatestMap(
    GridMapType& latestMap,
    const RobotPoseDeque& latestPoses,
    const ScanDataDeque& latestScanData)
{
    /* Vector for storing missed grid cell indices
     * Specified as static variable to reduce the performance loss
     * caused by memory allocations and deallocations */
    static std::vector<Point2D<int>> missedCellIndices;

    XAssert(!latestScanData.empty(),
            "Latest scan data should not be empty");
    XAssert(latestPoses.size() == latestScanData.size(),
            "The size of latest poses and latest scan data should be equal");

    /* Compute the scan points and bounding box */
    Point2D<double> minPos { std::numeric_limits<double>::max(),
                             std::numeric_limits<double>::max() };
    Point2D<double> maxPos { std::numeric_limits<double>::min(),
                             std::numeric_limits<double>::min() };
    std::vector<std::vector<Point2D<double>>> hitPoints;

    for (std::size_t i = 0; i < latestScanData.size(); ++i) {
        /* Retrieve the particle pose and scan data */
        const RobotPose2D<double>& particlePose = latestPoses[i];
        const auto& scanData = latestScanData[i];

        /* Compute the sensor pose from the particle pose */
        const RobotPose2D<double> sensorPose =
            Compound(particlePose, scanData->RelativeSensorPose());

        /* The latest grid map should contain the position of the sensor */
        minPos.mX = std::min(minPos.mX, sensorPose.mX);
        minPos.mY = std::min(minPos.mY, sensorPose.mY);
        maxPos.mX = std::max(maxPos.mX, sensorPose.mX);
        maxPos.mY = std::max(maxPos.mY, sensorPose.mY);

        /* Compute the minimum and maximum range of the scan */
        const double minRange = std::max(
            this->mMinUsableRange, scanData->MinRange());
        const double maxRange = std::min(
            this->mMaxUsableRange, scanData->MaxRange());

        /* Compute the scan points and bounding box */
        const std::size_t numOfScans = scanData->NumOfScans();
        std::vector<Point2D<double>> scanHitPoints;
        scanHitPoints.reserve(numOfScans);

        for (std::size_t j = 0; j < numOfScans; ++j) {
            const double range = scanData->RangeAt(j);

            if (range >= maxRange || range <= minRange)
                continue;

            /* Compute the hit point */
            Point2D<double> hitPoint = scanData->HitPoint(sensorPose, j);

            /* Update the bounding box */
            minPos.mX = std::min(minPos.mX, hitPoint.mX);
            minPos.mY = std::min(minPos.mY, hitPoint.mY);
            maxPos.mX = std::max(maxPos.mX, hitPoint.mX);
            maxPos.mY = std::max(maxPos.mY, hitPoint.mY);

            /* Append the hit point */
            scanHitPoints.emplace_back(std::move(hitPoint));
        }

        /* Push the hit points collected from the scan data */
        hitPoints.push_back(std::move(scanHitPoints));
    }

    /* Create a new grid map that contains all the scan points */
    latestMap.Expand(minPos.mX, maxPos.mX, minPos.mY, maxPos.mY);
    latestMap.Reset();

    /* Integrate the scan points into the grid map */
    for (std::size_t i = 0; i < latestScanData.size(); ++i) {
        /* Retrieve the particle pose and scan data */
        const RobotPose2D<double>& particlePose = latestPoses[i];
        const auto& scanData = latestScanData[i];

        /* Compute the sensor pose from the particle pose */
        const RobotPose2D<double> sensorPose =
            Compound(particlePose, scanData->RelativeSensorPose());
        /* Compute the grid cell index corresponding to the sensor pose */
        const Point2D<int> sensorCellIdx =
            latestMap.MapCoordinateToCellIndex(sensorPose.mX, sensorPose.mY);

        /* Integrate the scan points into the grid map */
        const std::size_t numOfHitPoints = hitPoints[i].size();

        for (std::size_t j = 0; j < numOfHitPoints; ++j) {
            /* Compute the grid cell index corresponding to the hit point */
            const Point2D<int> hitCellIdx =
                latestMap.MapCoordinateToCellIndex(hitPoints[i][j]);

            /* Compute the indices of the missed grid cells */
            this->ComputeMissedCellIndices(
                sensorCellIdx, hitCellIdx, missedCellIndices);

            /* Update the missed grid cells */
            const std::size_t numOfMissedCells = missedCellIndices.size();

            for (std::size_t k = 0; k < numOfMissedCells; ++k)
                latestMap.Update(missedCellIndices[k], this->mProbMiss);

            /* Update the hit grid cell */
            latestMap.Update(hitCellIdx, this->mProbHit);
        }
    }
}

/* Compute the indices of the missed cells using Bresenham algorithm */
void MapBuilder::ComputeMissedCellIndices(
    const Point2D<int>& startCellIdx,
    const Point2D<int>& endCellIdx,
    std::vector<Point2D<int>>& gridCellIndices) const
{
    /* Clear the grid cell indices */
    gridCellIndices.clear();
    /* Use Bresenham algorithm for computing indices */
    Bresenham(startCellIdx, endCellIdx, gridCellIndices);
    /* Remove the last item since it corresponds to the hit cell index */
    gridCellIndices.pop_back();
}

} /* namespace Mapping */
} /* namespace MyGMapping */
