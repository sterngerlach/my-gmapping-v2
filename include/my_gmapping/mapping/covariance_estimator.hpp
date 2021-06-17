
/* covariance_estimator.hpp */

#ifndef MY_GMAPPING_MAPPING_COVARIANCE_ESTIMATOR_HPP
#define MY_GMAPPING_MAPPING_COVARIANCE_ESTIMATOR_HPP

#include "my_gmapping/mapping/grid_map_types.hpp"
#include "my_gmapping/sensor/sensor_data.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>

namespace MyGMapping {
namespace Mapping {

/* Type definitions for convenience */
class CovarianceEstimator;
using CovarianceEstimatorPtr = std::unique_ptr<CovarianceEstimator>;

class CovarianceEstimator final
{
public:
    /* MapValues struct stores the four occupancy probability values
     * at integer coordinates closest to the specified grid cell indices
     * in floating-point, which are later used in bilinear interpolations
     * and gradient (with respect to the robot pose) computations */
    struct MapValues
    {
        /* Constructor */
        MapValues(const double dx, const double dy,
                  const double m00, const double m01,
                  const double m10, const double m11);

        /* Compute the smoothed map value using bilinear interpolation */
        double BilinearInterpolation() const;

        /* Relative position of the specified indices */
        const double mDeltaX;
        const double mDeltaY;
        /* Four neighboring map values */
        const double mValue00;
        const double mValue01;
        const double mValue10;
        const double mValue11;
    };

public:
    /* Constructor */
    CovarianceEstimator(const double covarianceScale) :
        mCovarianceScale(covarianceScale) { }
    /* Destructor */
    ~CovarianceEstimator() = default;

    /* Compute a covariance matrix in a world coordinate frame */
    Eigen::Matrix3d ComputeCovariance(
        const GridMapInterface& gridMap,
        const Sensor::ScanDataPtr<double>& scanData,
        const RobotPose2D<double>& sensorPose) const;

    /* Compute a Hessian matrix using a gradient of a smoothed map function
     * with respect to the sensor pose in a world coordinate frame */
    Eigen::Matrix3d ComputeHessian(
        const GridMapInterface& gridMap,
        const Sensor::ScanDataPtr<double>& scanData,
        const RobotPose2D<double>& sensorPose) const;

    /* Compute a gradient of the smoothed map function
     * at the specified scan point with respect to the robot pose */
    Eigen::Vector3d ComputeScaledMapGradSensorPose(
        const MapValues mapValues,
        const Point2D<double>& rotatedScanPoint) const;

    /* Compute a (normalized) gradient of the smoothed map function
     * based on the bilinear interpolation with respect to the
     * specified position of the map */
    Eigen::Vector2d ComputeScaledMapGradMapPoint(
        const MapValues mapValues) const;

    /* Get the four occupancy probability values at the integer coordinates
     * closest to the specified grid cell indices in floating-point */
    MapValues GetClosestMapValues(
        const GridMapInterface& gridMap,
        const Point2D<double>& gridCellIdx) const;

private:
    /* Scale factor for computing the covariance matrix */
    const double mCovarianceScale;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_MAPPING_COVARIANCE_ESTIMATOR_HPP */
