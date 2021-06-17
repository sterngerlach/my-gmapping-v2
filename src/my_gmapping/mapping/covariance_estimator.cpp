
/* covariance_estimator.cpp */

#include "my_gmapping/mapping/covariance_estimator.hpp"

namespace MyGMapping {
namespace Mapping {

/* Constructor */
CovarianceEstimator::MapValues::MapValues(
    const double dx, const double dy,
    const double m00, const double m01,
    const double m10, const double m11) :
    mDeltaX(dx), mDeltaY(dy),
    mValue00(m00), mValue01(m01),
    mValue10(m10), mValue11(m11)
{
    /* Check the relative position */
    Assert(this->mDeltaX >= 0.0 && this->mDeltaX <= 1.0);
    Assert(this->mDeltaY >= 0.0 && this->mDeltaY <= 1.0);
}

/* Compute the smoothed map value using bilinear interpolation */
double CovarianceEstimator::MapValues::BilinearInterpolation() const
{
    const double interpolated =
        this->mDeltaY * (this->mDeltaX * this->mValue11 +
                         (1.0 - this->mDeltaX) * this->mValue01) +
        (1.0 - this->mDeltaY) * (this->mDeltaX * this->mValue10 +
                                 (1.0 - this->mDeltaX) * this->mValue00);

    return interpolated;
}

/* Compute a covariance matrix in a world coordinate frame */
Eigen::Matrix3d CovarianceEstimator::ComputeCovariance(
    const GridMapInterface& gridMap,
    const Sensor::ScanDataPtr<double>& scanData,
    const RobotPose2D<double>& sensorPose) const
{
    /* Compute the Hessian matrix */
    const Eigen::Matrix3d hessianMat =
        this->ComputeHessian(gridMap, scanData, sensorPose);
    /* Compute the covariance matrix by scaling the inverse of Hessian */
    const Eigen::Matrix3d covarianceMat =
        hessianMat.inverse() * this->mCovarianceScale;

    return covarianceMat;
}

/* Compute a Hessian matrix using a gradient of a smoothed map function
 * with respect to the sensor pose in a world coordinate frame */
Eigen::Matrix3d CovarianceEstimator::ComputeHessian(
    const GridMapInterface& gridMap,
    const Sensor::ScanDataPtr<double>& scanData,
    const RobotPose2D<double>& sensorPose) const
{
    /* Compute an approximated Hessian matrix and a residual vector */
    Eigen::Matrix3d hessianMat = Eigen::Matrix3d::Zero();

    const double reciprocalResolution = 1.0 / gridMap.Resolution();
    const std::size_t numOfScans = scanData->NumOfScans();

    for (std::size_t i = 0; i < numOfScans; ++i) {
        /* Compute the grid cell index corresponding to the scan point */
        const Point2D<double> hitPoint =
            scanData->HitPoint(sensorPose, i);
        const Point2D<double> floatingIdx =
            gridMap.PositionToIndexF(hitPoint.mX, hitPoint.mY);
        /* Get the map values closest to the floating-point index */
        const auto mapValues = this->GetClosestMapValues(gridMap, floatingIdx);

        /* Evaluate the (scaled) gradient of the smoothed map function
         * at the scan point with respect to the sensor pose */
        const Point2D<double> rotatedScanPoint {
            hitPoint.mX - sensorPose.mX,
            hitPoint.mY - sensorPose.mY };
        const Eigen::Vector3d scaledMapGrad =
            this->ComputeScaledMapGradSensorPose(mapValues, rotatedScanPoint);
        const Eigen::Vector3d mapGrad = scaledMapGrad * reciprocalResolution;

        /* Compute the Hessian matrix for the scan point */
        const Eigen::Matrix3d scanHessianMat = mapGrad * mapGrad.transpose();
        hessianMat += scanHessianMat;
    }

    return hessianMat;
}

/* Compute a gradient of the smoothed map function
 * at the specified scan point with respect to the robot pose */
Eigen::Vector3d CovarianceEstimator::ComputeScaledMapGradSensorPose(
    const MapValues mapValues,
    const Point2D<double>& rotatedScanPoint) const
{
    /* Compute a gradient of the smoothed occupancy probability value
     * with respect to the scan point */
    const Eigen::Vector2d scaledMapGrad =
        this->ComputeScaledMapGradMapPoint(mapValues);

    /* Compute a scaled gradient of the smoothed occupancy probability value
     * at the scan point with respect to the robot pose */
    /* Divide by the grid map resolution to the following values
     * to get the actual gradient */
    const double gradX = scaledMapGrad(0);
    const double gradY = scaledMapGrad(1);
    const double gradTheta = -rotatedScanPoint.mY * scaledMapGrad(0) +
                              rotatedScanPoint.mX * scaledMapGrad(1);

    return Eigen::Vector3d { gradX, gradY, gradTheta };
}

/* Compute a (normalized) gradient of the smoothed map function
 * based on the bilinear interpolation with respect to the
 * specified map-local position of the map */
Eigen::Vector2d CovarianceEstimator::ComputeScaledMapGradMapPoint(
    const MapValues mapValues) const
{
    /* Refer to the Equation (5) and (6) in the following paper:
     * Stefan Kohlbrecher, Oskar von Stryk, Johannes Meyer and Uwe Klingauf.
     * "A Flexible and Scalable SLAM System with Full 3D Motion Estimation,"
     * in the Proceedings of the IEEE International Symposium on Safety,
     * Security and Rescue Robotics (SSRR), 2011. */

    /* Divide by the grid map resolution to the following values
     * to get the actual gradient */
    const double scaledGradX =
        mapValues.mDeltaY * (mapValues.mValue11 - mapValues.mValue01) +
        (1.0 - mapValues.mDeltaY) * (mapValues.mValue10 - mapValues.mValue00);
    const double scaledGradY =
        mapValues.mDeltaX * (mapValues.mValue11 - mapValues.mValue10) +
        (1.0 - mapValues.mDeltaX) * (mapValues.mValue01 - mapValues.mValue00);

    return Eigen::Vector2d { scaledGradX, scaledGradY };
}

/* Get the four occupancy probability values at the integer coordinates
 * closest to the specified grid cell indices in floating-point */
CovarianceEstimator::MapValues CovarianceEstimator::GetClosestMapValues(
    const GridMapInterface& gridMap,
    const Point2D<double>& gridCellIdx) const
{
    /* Obtain the closest integer coordinates */
    const double x0 = std::floor(gridCellIdx.mX);
    const double y0 = std::floor(gridCellIdx.mY);
    const double dx = gridCellIdx.mX - x0;
    const double dy = gridCellIdx.mY - y0;

    /* Clamp the integer coordinates, since the grid cell index
     * could be out-of-bounds */
    const int xc0 = std::max(static_cast<int>(x0), 0);
    const int yc0 = std::max(static_cast<int>(y0), 0);
    const int xc1 = std::min(xc0 + 1, gridMap.Cols() - 1);
    const int yc1 = std::min(yc0 + 1, gridMap.Rows() - 1);

    /* Obtain the occupancy probability values at four integer coordinates */
    const double m00 = gridMap.ProbabilityOr(yc0, xc0, 0.5);
    const double m01 = gridMap.ProbabilityOr(yc1, xc0, 0.5);
    const double m10 = gridMap.ProbabilityOr(yc0, xc1, 0.5);
    const double m11 = gridMap.ProbabilityOr(yc1, xc1, 0.5);

    return MapValues { dx, dy, m00, m01, m10, m11 };
}

} /* namespace Mapping */
} /* namespace MyGMapping */
