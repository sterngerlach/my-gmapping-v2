
/* sensor_data.hpp */

#ifndef MY_GMAPPING_SENSOR_DATA_HPP
#define MY_GMAPPING_SENSOR_DATA_HPP

#include <memory>
#include <vector>

#include "my_gmapping/pose.hpp"

namespace MyGMapping {
namespace Sensor {

class SensorData
{
public:
    /* Constructor */
    SensorData(const std::string& sensorId, double timeStamp) :
        mSensorId(sensorId), mTimeStamp(timeStamp) { }
    /* Destructor */
    virtual ~SensorData() = default;

    /* Retrieve sensor id string */
    inline const std::string& SensorId() const { return this->mSensorId; }
    /* Retrieve timestamp */
    inline double TimeStamp() const { return this->mTimeStamp; }

protected:
    std::string     mSensorId;
    double          mTimeStamp;
};

template <typename T>
class OdometryData final : public SensorData
{
public:
    /* Constructor */
    OdometryData(const std::string& sensorId,
                 double timeStamp,
                 const RobotPose2D<T>& pose,
                 const RobotPose2D<T>& velocity) :
        SensorData(sensorId, timeStamp),
        mPose(pose),
        mVelocity(velocity) { }

    /* Destructor */
    ~OdometryData() = default;

    /* Retrieve the odometry pose in base frame */
    inline const RobotPose2D<T>& Pose() const { return this->mPose; }
    /* Retrieve the velocity in base frame */
    inline const RobotPose2D<T>& Velocity() const { return this->mVelocity; }

private:
    /* Odometry pose in base frame */
    RobotPose2D<T> mPose;
    /* Velocity in base frame */
    RobotPose2D<T> mVelocity;
};

template <typename T>
class ScanData final : public SensorData
{
public:
    /* Constructor */
    ScanData(const std::string& sensorId,
             double timeStamp,
             const RobotPose2D<T>& odomPose,
             const RobotPose2D<T>& velocity,
             const RobotPose2D<T>& relativeSensorPose,
             const T& minRange,
             const T& maxRange,
             const T& minAngle,
             const T& maxAngle,
             std::vector<T>&& angles,
             std::vector<T>&& ranges) :
        SensorData(sensorId, timeStamp),
        mOdomPose(odomPose),
        mVelocity(velocity),
        mRelativeSensorPose(relativeSensorPose),
        mMinRange(minRange),
        mMaxRange(maxRange),
        mMinAngle(minAngle),
        mMaxAngle(maxAngle),
        mAngles(angles),
        mRanges(ranges) { }
    
    /* Destructor */
    ~ScanData() = default;

    /* Retrieve the odometry pose in base frame */
    inline const RobotPose2D<T>& OdomPose() const
    { return this->mOdomPose; }
    /* Retrieve the robot velocity in base frame */
    inline const RobotPose2D<T>& Velocity() const
    { return this->mVelocity; }

    /* Retrieve the sensor pose in robot frame */
    inline const RobotPose2D<T>& RelativeSensorPose() const
    { return this->mRelativeSensorPose; }

    /* Retrieve the minimum range in meters */
    inline T MinRange() const { return this->mMinRange; }
    /* Retrieve the maximum range in meters */
    inline T MaxRange() const { return this->mMaxRange; }
    /* Retrieve the minimum angle in radians */
    inline T MinAngle() const { return this->mMinAngle; }
    /* Retrieve the maximum angle in radians */
    inline T MaxAngle() const { return this->mMaxAngle; }

    /* Retrieve the angle data */
    inline const std::vector<T>& Angles() const { return this->mAngles; }
    /* Retrieve the range data */
    inline const std::vector<T>& Ranges() const { return this->mRanges; }
    /* Retrieve the number of scans */
    inline std::size_t NumOfScans() const { return this->mRanges.size(); }

    /* Retrieve the scan range at the specified index */
    inline T RangeAt(const std::size_t scanIdx) const
    { return this->mRanges.at(scanIdx); }
    /* Retrieve the scan angle at the specified index */
    inline T AngleAt(const std::size_t scanIdx) const
    { return this->mAngles.at(scanIdx); }

    /* Compute the hit point of the specified scan */
    Point2D<double> HitPoint(
        const RobotPose2D<double>& sensorPose,
        const std::size_t scanIdx) const;
    /* Compute the hit pose of the specified scan */
    RobotPose2D<double> HitPose(
        const RobotPose2D<double>& sensorPose,
        const std::size_t scanIdx) const;

    /* Compute the hit point of the specified scan in a sensor-local frame */
    Point2D<double> SensorLocalHitPoint(
        const std::size_t scanIdx) const;
    /* Compute the hit pose of the specified scan in a sensor-local frame */
    RobotPose2D<double> SensorLocalHitPose(
        const std::size_t scanIdx) const;

    /* Compute the hit point and missed point of the specified scan */
    void HitAndMissedPoint(
        const RobotPose2D<double>& sensorPose,
        const std::size_t scanIdx,
        const T hitAndMissedDist,
        Point2D<double>& hitPoint,
        Point2D<double>& missedPoint) const;

private:
    /* Odometry pose in base frame (required) */
    RobotPose2D<T>  mOdomPose;
    /* Robot velocity in base frame (optional) */
    RobotPose2D<T>  mVelocity;
    /* Relative sensor pose in robot frame */
    RobotPose2D<T>  mRelativeSensorPose;
    /* Minimum range in meters */
    T               mMinRange;
    /* Maximum range in meters */
    T               mMaxRange;
    /* Minimum angle in radians */
    T               mMinAngle;
    /* Maximum angle in radians */
    T               mMaxAngle;
    /* Angle data */
    std::vector<T>  mAngles;
    /* Range data */
    std::vector<T>  mRanges;
};

/* Compute the hit point of the specified scan */
template <typename T>
Point2D<double> ScanData<T>::HitPoint(
    const RobotPose2D<double>& sensorPose,
    const std::size_t scanIdx) const
{
    const T scanRange = this->RangeAt(scanIdx);
    const T scanAngle = this->AngleAt(scanIdx);
    const T cosTheta = std::cos(sensorPose.mTheta + scanAngle);
    const T sinTheta = std::sin(sensorPose.mTheta + scanAngle);

    return Point2D<double> {
        sensorPose.mX + scanRange * cosTheta,
        sensorPose.mY + scanRange * sinTheta };
}

/* Compute the hit pose of the specified scan */
template <typename T>
RobotPose2D<double> ScanData<T>::HitPose(
    const RobotPose2D<double>& sensorPose,
    const std::size_t scanIdx) const
{
    const T scanRange = this->RangeAt(scanIdx);
    const T scanAngle = this->AngleAt(scanIdx);
    const T cosTheta = std::cos(sensorPose.mTheta + scanAngle);
    const T sinTheta = std::sin(sensorPose.mTheta + scanAngle);

    return RobotPose2D<double> {
        sensorPose.mX + scanRange * cosTheta,
        sensorPose.mY + scanRange * sinTheta,
        sensorPose.mTheta + scanAngle };
}

/* Compute the hit point of the specified scan in a sensor-local frame */
template <typename T>
Point2D<double> ScanData<T>::SensorLocalHitPoint(
    const std::size_t scanIdx) const
{
    const T scanRange = this->RangeAt(scanIdx);
    const T scanAngle = this->AngleAt(scanIdx);
    const T cosTheta = std::cos(scanAngle);
    const T sinTheta = std::sin(scanAngle);

    return Point2D<double> {
        scanRange * cosTheta, scanRange * sinTheta };
}

/* Compute the hit pose of the specified scan in a sensor-local frame */
template <typename T>
RobotPose2D<double> ScanData<T>::SensorLocalHitPose(
    const std::size_t scanIdx) const
{
    const T scanRange = this->RangeAt(scanIdx);
    const T scanAngle = this->AngleAt(scanIdx);
    const T cosTheta = std::cos(scanAngle);
    const T sinTheta = std::sin(scanAngle);

    return RobotPose2D<double> {
        scanRange * cosTheta, scanAngle * sinTheta, scanAngle };
}

/* Compute the hit point and missed point of the specified scan */
template <typename T>
void ScanData<T>::HitAndMissedPoint(
    const RobotPose2D<double>& sensorPose,
    const std::size_t scanIdx,
    const T hitAndMissedDist,
    Point2D<double>& hitPoint,
    Point2D<double>& missedPoint) const
{
    const T scanRange = this->RangeAt(scanIdx);
    const T scanAngle = this->AngleAt(scanIdx);
    const T cosTheta = std::cos(sensorPose.mTheta + scanAngle);
    const T sinTheta = std::sin(sensorPose.mTheta + scanAngle);

    /* Compute the hit point */
    hitPoint.mX = sensorPose.mX + scanRange * cosTheta;
    hitPoint.mY = sensorPose.mY + scanRange * sinTheta;

    /* Compute the missed point */
    missedPoint.mX = sensorPose.mX + (scanRange - hitAndMissedDist) * cosTheta;
    missedPoint.mY = sensorPose.mY + (scanRange - hitAndMissedDist) * sinTheta;

    return;
}

/* Type definitions */
using SensorDataPtr = std::shared_ptr<SensorData>;

template <typename T>
using OdometryDataPtr = std::shared_ptr<OdometryData<T>>;

template <typename T>
using ScanDataPtr = std::shared_ptr<ScanData<T>>;

} /* namespace Sensor */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_SENSOR_DATA_HPP */
