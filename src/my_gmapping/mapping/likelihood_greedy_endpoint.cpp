
/* likelihood_greedy_endpoint.cpp */

#include "my_gmapping/mapping/likelihood_greedy_endpoint.hpp"

namespace MyGMapping {
namespace Mapping {

/* Constructor */
LikelihoodGreedyEndpoint::LikelihoodGreedyEndpoint(
    const double mapResolution,
    const double minUsableRange,
    const double maxUsableRange,
    const double hitAndMissedCellDist,
    const double occupancyThreshold,
    const double gaussianSigma,
    const int kernelSize,
    const double likelihoodScale,
    const bool useLocalMap,
    const int localMapSize) :
    LikelihoodFunction(),
    mMapResolution(mapResolution),
    mMinUsableRange(minUsableRange),
    mMaxUsableRange(maxUsableRange),
    mHitAndMissedCellDist(hitAndMissedCellDist),
    mOccupancyThreshold(occupancyThreshold),
    mGaussianSigma(gaussianSigma),
    mKernelSize(kernelSize),
    mLikelihoodScale(likelihoodScale),
    mLikelihoodTable(nullptr),
    mDefaultLikelihood(0.0),
    mUseLocalMap(useLocalMap),
    mLocalMapSize(localMapSize)
{
    /* Compute the lookup table for the likelihood value */
    this->SetupLookupTable();
}

/* Calculate observation likelihood based on squared distance
 * between scan point and matched cell position */
double LikelihoodGreedyEndpoint::Likelihood(
    const GridMapInterface& gridMap,
    const Sensor::ScanDataPtr<double>& scanData,
    const RobotPose2D<double>& sensorPose)
{
    double likelihoodSum = 0.0;

    const double minRange = std::max(
        this->mMinUsableRange, scanData->MinRange());
    const double maxRange = std::min(
        this->mMaxUsableRange, scanData->MaxRange());

    /* Compute the grid cell index corresponding to the sensor pose */
    const Point2D<int> sensorPoseIdx =
        gridMap.PositionToIndex(sensorPose.mX, sensorPose.mY);

    /* Compute the grid cell index range of the local map */
    const int localMapRadius = this->mLocalMapSize / 2;
    const Point2D<int> localMapIdxMin {
        sensorPoseIdx.mX - localMapRadius, sensorPoseIdx.mY - localMapRadius };
    const Point2D<int> localMapIdxMax {
        sensorPoseIdx.mX + localMapRadius, sensorPoseIdx.mY + localMapRadius };

    for (std::size_t i = 0; i < scanData->NumOfScans(); ++i) {
        const double range = scanData->RangeAt(i);

        if (range >= maxRange || range <= minRange)
            continue;

        /* Compute the grid cell indices corresponding to the scan point and
         * the missed point */
        Point2D<double> hitPoint;
        Point2D<double> missedPoint;
        scanData->HitAndMissedPoint(sensorPose, i, this->mHitAndMissedCellDist,
                                    hitPoint, missedPoint);

        const Point2D<int> hitIdx =
            gridMap.PositionToIndex(hitPoint.mX, hitPoint.mY);
        const Point2D<int> missedIdx =
            gridMap.PositionToIndex(missedPoint.mX, missedPoint.mY);

        /* Search the best cell index from the window */
        const double unknownProb = gridMap.UnknownProbability();
        double likelihoodMax = this->mDefaultLikelihood;

        for (int ky = -this->mKernelSize; ky <= this->mKernelSize; ++ky) {
            for (int kx = -this->mKernelSize; kx <= this->mKernelSize; ++kx) {
                const double hitProb = gridMap.ProbabilityOr(
                    hitIdx.mY + ky, hitIdx.mX + kx, unknownProb);
                const double missProb = gridMap.ProbabilityOr(
                    missedIdx.mY + ky, missedIdx.mX + kx, unknownProb);

                /* Check if the cell is inside of the local map */
                const bool hitInside = hitIdx.mX >= localMapIdxMin.mX &&
                                       hitIdx.mY >= localMapIdxMin.mY &&
                                       hitIdx.mX < localMapIdxMax.mX &&
                                       hitIdx.mY < localMapIdxMax.mY;
                const bool missedInside = missedIdx.mX >= localMapIdxMin.mX &&
                                          missedIdx.mY >= localMapIdxMin.mY &&
                                          missedIdx.mX < localMapIdxMax.mX &&
                                          missedIdx.mY < localMapIdxMax.mY;

                /* Skip if the cell index is out of bounds */
                if (this->mUseLocalMap && (!hitInside || !missedInside))
                    continue;

                /* Skip if cell contains unknown occupancy probability */
                if (!this->mUseLocalMap &&
                    (hitProb == unknownProb || missProb == unknownProb))
                    continue;

                /* Skip if the occupancy probability of the cell
                 * that is assumed to be hit is less than the threshold or
                 * the occupancy probability of the missed cell is greater
                 * than the threshold */
                if (hitProb < this->mOccupancyThreshold ||
                    missProb > this->mOccupancyThreshold)
                    continue;

                /* Retrieve the likelihood value using the lookup table */
                const int idxX = this->mKernelSize + kx;
                const int idxY = this->mKernelSize + ky;
                const int tableIdx = idxY * (2 * this->mKernelSize + 1) + idxX;
                const double likelihood = this->mLikelihoodTable[tableIdx];

                /* Update the maximum likelihood value */
                likelihoodMax = std::max(likelihoodMax, likelihood);
            }
        }

        /* Add to the scan matching score value
         * Score represents the log-likelihood of the observation */
        likelihoodSum += likelihoodMax;
    }

    /* Scale the observation likelihood
     * to prevent underflow in weight computation */
    likelihoodSum *= this->mLikelihoodScale;

    return likelihoodSum;
}

/* Setup the lookup table for the likelihood value */
void LikelihoodGreedyEndpoint::SetupLookupTable()
{
    /* Allocate the lookup table */
    const int kernelSize = 2 * this->mKernelSize + 1;
    this->mLikelihoodTable.reset(new double[kernelSize * kernelSize]);

    /* Compute the likelihood values in the lookup table */
    for (int ky = -this->mKernelSize; ky <= this->mKernelSize; ++ky) {
        for (int kx = -this->mKernelSize; kx <= this->mKernelSize; ++kx) {
            /* Compute the squared distance using the offsets `kx` and `ky` */
            const double diffX = this->mMapResolution * kx;
            const double diffY = this->mMapResolution * ky;
            const double squaredDist = diffX * diffX + diffY * diffY;

            /* Compute the negative log-likelihood of the observation */
            const double likelihoodValue =
                -squaredDist / (0.5 * this->mGaussianSigma);

            /* Set the likelihood value to the lookup table */
            const int idxX = this->mKernelSize + kx;
            const int idxY = this->mKernelSize + ky;
            const int tableIdx = idxY * kernelSize + idxX;
            this->mLikelihoodTable[tableIdx] = likelihoodValue;
        }
    }

    /* Compute the default likelihood value */
    const double maxDiffX = this->mMapResolution * (this->mKernelSize + 1);
    const double maxDiffY = this->mMapResolution * (this->mKernelSize + 1);
    const double maxSquaredDist = maxDiffX * maxDiffX + maxDiffY * maxDiffY;

    /* Compute the negative log-likelihood of the observation */
    this->mDefaultLikelihood = -maxSquaredDist / (0.5 * this->mGaussianSigma);
}

} /* namespace Mapping */
} /* namespace MyGMapping */
