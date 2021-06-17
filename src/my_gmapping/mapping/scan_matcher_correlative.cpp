
/* scan_matcher_correlative.cpp */

#include "my_gmapping/mapping/scan_matcher_correlative.hpp"
#include "my_gmapping/mapping/grid_map_builder.hpp"

namespace MyGMapping {
namespace Mapping {

/* Constructor */
ScanMatcherCorrelativeMetrics::ScanMatcherCorrelativeMetrics(
    const std::string& scanMatcherName) :
    mInputSetupTime(nullptr),
    mOptimizationTime(nullptr),
    mDiffTranslation(nullptr),
    mDiffRotation(nullptr),
    mWinSizeX(nullptr),
    mWinSizeY(nullptr),
    mWinSizeTheta(nullptr),
    mStepSizeX(nullptr),
    mStepSizeY(nullptr),
    mStepSizeTheta(nullptr),
    mNumOfIgnoredNodes(nullptr),
    mNumOfProcessedNodes(nullptr),
    mScoreValue(nullptr),
    mLikelihoodValue(nullptr),
    mNumOfScans(nullptr)
{
    /* Retrieve the metrics manager instance */
    auto* const pMetricManager = Metric::MetricManager::Instance();

    /* Register the counter metrics */
    this->mNumOfIgnoredNodes = pMetricManager->AddCounter(
        scanMatcherName + ".NumOfIgnoredNodes");
    this->mNumOfProcessedNodes = pMetricManager->AddCounter(
        scanMatcherName + ".NumOfProcessedNodes");

    /* Register the distribution metrics */
    this->mInputSetupTime = pMetricManager->AddDistribution(
        scanMatcherName + ".InputSetupTime");
    this->mOptimizationTime = pMetricManager->AddDistribution(
        scanMatcherName + ".OptimizationTime");
    this->mDiffTranslation = pMetricManager->AddDistribution(
        scanMatcherName + ".DiffTranslation");
    this->mDiffRotation = pMetricManager->AddDistribution(
        scanMatcherName + ".DiffRotation");

    this->mWinSizeX = pMetricManager->AddDistribution(
        scanMatcherName + ".WinSizeX");
    this->mWinSizeY = pMetricManager->AddDistribution(
        scanMatcherName + ".WinSizeY");
    this->mWinSizeTheta = pMetricManager->AddDistribution(
        scanMatcherName + ".WinSizeTheta");
    this->mStepSizeX = pMetricManager->AddDistribution(
        scanMatcherName + ".StepSizeX");
    this->mStepSizeY = pMetricManager->AddDistribution(
        scanMatcherName + ".StepSizeY");
    this->mStepSizeTheta = pMetricManager->AddDistribution(
        scanMatcherName + ".StepSizeTheta");

    /* Register the value sequence metrics */
    this->mScoreValue = pMetricManager->AddValueSequence<float>(
        scanMatcherName + ".ScoreValue");
    this->mLikelihoodValue = pMetricManager->AddValueSequence<float>(
        scanMatcherName + ".LikelihoodValue");
    this->mNumOfScans = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".NumOfScans");
}

/* Constructor */
ScanMatcherCorrelative::ScanMatcherCorrelative(
    const std::string& scanMatcherName,
    LikelihoodFuncPtr&& likelihoodFunc,
    const bool useCroppedMap,
    const int croppedMapSizeX,
    const int croppedMapSizeY,
    const int lowResolution,
    const double rangeX,
    const double rangeY,
    const double rangeTheta,
    const double scanRangeMax) :
    ScanMatcher(scanMatcherName),
    mLikelihoodFunc(std::move(likelihoodFunc)),
    mUseCroppedMap(useCroppedMap),
    mCroppedMapSizeX(croppedMapSizeX),
    mCroppedMapSizeY(croppedMapSizeY),
    mLowResolution(lowResolution),
    mRangeX(rangeX),
    mRangeY(rangeY),
    mRangeTheta(rangeTheta),
    mScanRangeMax(scanRangeMax),
    mMetrics(scanMatcherName)
{
}

/* Optimize the particle poses based on the correlative scan matching */
void ScanMatcherCorrelative::OptimizePose(
    const std::size_t numOfParticles,
    const std::vector<const GridMap*>& particleMaps,
    const Sensor::ScanDataPtr<double>& scanData,
    const std::vector<RobotPose2D<double>>& initialPoses,
    std::vector<RobotPose2D<double>>& estimatedPoses,
    std::vector<double>& likelihoodValues)
{
    /* Input checks */
    assert(particleMaps.size() == numOfParticles);
    assert(initialPoses.size() == numOfParticles);

    estimatedPoses.clear();
    estimatedPoses.resize(numOfParticles);
    likelihoodValues.clear();
    likelihoodValues.resize(numOfParticles);

    std::vector<double> normalizedLikelihoods;
    normalizedLikelihoods.resize(numOfParticles);
    std::vector<double> normalizedScores;
    normalizedScores.resize(numOfParticles);

    /* Optimize the pose for each particle (parallelized using OpenMP) */
#pragma omp parallel for
    for (std::size_t i = 0; i < numOfParticles; ++i) {
        /* Create the timer for each thread */
        Metric::Timer timer;

        /* Create the grid map to store the intermediate result and the final
         * result (low-resolution coarse grid map) */
        ConstMap intermediateMap {
            particleMaps[i]->Resolution(), particleMaps[i]->BlockSize(),
            this->mCroppedMapSizeY, this->mCroppedMapSizeX };
        ConstMap precompMap {
            particleMaps[i]->Resolution(), particleMaps[i]->BlockSize(),
            this->mCroppedMapSizeY, this->mCroppedMapSizeX };

        /* Determine the desired size of the grid map for scan matching */
        const int mapColsMax = this->mUseCroppedMap ?
            this->mCroppedMapSizeX : particleMaps[i]->Cols();
        const int mapRowsMax = this->mUseCroppedMap ?
            this->mCroppedMapSizeY : particleMaps[i]->Rows();

        /* Compute the center position of the cropped grid map */
        const Point2D<int> desiredCenterIdx = particleMaps[i]->PositionToIndex(
            initialPoses[i].mX, initialPoses[i].mY);
        /* Use std::clamp() here since the above `centerIdx` could be
         * out-of-bounds (the current particle is outside of the grid map) */
        const Point2D<int> possibleIdxMin {
            std::clamp(desiredCenterIdx.mX - mapColsMax, 0,
                       particleMaps[i]->Cols() - 1),
            std::clamp(desiredCenterIdx.mY - mapRowsMax, 0,
                       particleMaps[i]->Rows() - 1) };
        const Point2D<int> possibleIdxMax {
            std::clamp(desiredCenterIdx.mX + mapColsMax, 1,
                       particleMaps[i]->Cols()),
            std::clamp(desiredCenterIdx.mY + mapRowsMax, 1,
                       particleMaps[i]->Rows()) };
        const Point2D<int> centerIdx {
            (possibleIdxMax.mX + possibleIdxMin.mX) / 2,
            (possibleIdxMax.mY + possibleIdxMin.mY) / 2 };

        /* Crop the grid map around the center position */
        const Point2D<int> idxMin {
            std::max(centerIdx.mX - mapColsMax / 2, 0),
            std::max(centerIdx.mY - mapRowsMax / 2, 0) };
        const Point2D<int> idxMax {
            std::min(centerIdx.mX + mapColsMax / 2, particleMaps[i]->Cols()),
            std::min(centerIdx.mY + mapRowsMax / 2, particleMaps[i]->Rows()) };
        const BoundingBox<int> boundingBox { idxMin, idxMax };

        /* Precompute the coarser grid map from the cropped grid map */
        PrecomputeGridMap(*particleMaps[i], this->mLowResolution, idxMin,
                          intermediateMap, precompMap);

        /* Update the processing time for setting up the input */
        this->mMetrics.mInputSetupTime->Observe(timer.ElapsedMicro());
        timer.Start();

        /* Optimize the particle pose based on the correlative scan matching */
        RobotPose2D<double> estimatedPose;
        double logLikelihood;
        double normalizedLikelihood;
        double normalizedScore;
        this->OptimizePoseCore(*particleMaps[i], boundingBox, precompMap,
                               scanData, initialPoses[i], 0.0,
                               estimatedPose, logLikelihood,
                               normalizedLikelihood, normalizedScore);

        /* Update the processing time for the optimization */
        this->mMetrics.mOptimizationTime->Observe(timer.ElapsedMicro());
        timer.Stop();

        /* Set the refined particle pose and the observation likelihood */
        estimatedPoses[i] = estimatedPose;
        likelihoodValues[i] = logLikelihood;
        normalizedLikelihoods[i] = normalizedLikelihood;
        normalizedScores[i] = normalizedScore;
    }

    /* Determine the maximum likelihood value and its corresponding score */
    const auto bestIt = std::max_element(
        likelihoodValues.begin(), likelihoodValues.end());
    const auto bestIdx = std::distance(likelihoodValues.begin(), bestIt);

    /* Update the normalized score value of the best particle */
    this->mMetrics.mScoreValue->Observe(normalizedScores[bestIdx]);
    /* Update the normalized likelihood value of the best particle */
    this->mMetrics.mLikelihoodValue->Observe(normalizedLikelihoods[bestIdx]);
    /* Update the number of the scan points */
    this->mMetrics.mNumOfScans->Observe(scanData->NumOfScans());
}

/* Optimize the particle pose based on the correlative scan matching */
void ScanMatcherCorrelative::OptimizePoseCore(
    const GridMapInterface& gridMap,
    const BoundingBox<int>& boundingBox,
    const GridMapInterface& coarseGridMap,
    const Sensor::ScanDataPtr<double>& scanData,
    const RobotPose2D<double>& initialPose,
    const double normalizedScoreThreshold,
    RobotPose2D<double>& estimatedPose,
    double& likelihoodValue,
    double& normalizedLikelihood,
    double& normalizedScore)
{
    /* Find the best particle pose from the search window */
    /* Compute the sensor pose from the initial particle pose */
    const RobotPose2D<double> sensorPose =
        Compound(initialPose, scanData->RelativeSensorPose());

    /* Determine the search step */
    double stepX;
    double stepY;
    double stepTheta;
    this->ComputeSearchStep(gridMap, scanData, stepX, stepY, stepTheta);

    /* Determine the search window */
    /* 'winX' and 'winY' are in the number of grid cells */
    const int winX = static_cast<int>(
        std::ceil(0.5 * this->mRangeX / stepX));
    const int winY = static_cast<int>(
        std::ceil(0.5 * this->mRangeY / stepY));
    const int winTheta = static_cast<int>(
        std::ceil(0.5 * this->mRangeTheta / stepTheta));

    /* Perform the scan matching against the low resolution grid map */
    double scoreMax = normalizedScoreThreshold;
    int bestWinX = -winX;
    int bestWinY = -winY;
    int bestWinTheta = -winTheta;
    /* Setup the number of the skipped nodes and the processed nodes */
    int numOfIgnoredNodes = 0;
    int numOfProcessedNodes = 0;

    /* Compute the grid cell indices for scan points */
    std::vector<Point2D<int>> scanIndices;
    scanIndices.reserve(scanData->NumOfScans());

    for (int t = -winTheta; t <= winTheta; ++t) {
        /* Compute the grid cell indices for scan points */
        const RobotPose2D<double> currentSensorPose {
            sensorPose.mX, sensorPose.mY, sensorPose.mTheta + stepTheta * t };
        this->ComputeScanIndices(
            gridMap, currentSensorPose, scanData, scanIndices);

        /* 'winX' and 'winY' are represented in the number of grid cells */
        /* For given 't', the projected scan points 'scanIndices' are
         * related by pure translation for the 'x' and 'y' search directions */
        for (int x = -winX; x < winX; x += this->mLowResolution) {
            for (int y = -winY; y < winY; y += this->mLowResolution) {
                /* Evaluate the normalized matching score
                 * Add index offsets `gridMapMinIdx.mX` and `gridMapMinIdx.mY`
                 * since the grid cell indices for the scan points `scanIndices`
                 * are computed from the original grid map `gridMap` */
                const BoundingBox<int>& coarseMapBox {
                    0, 0, coarseGridMap.Cols(), coarseGridMap.Rows() };
                const double normalizedScore = this->ComputeScore(
                    coarseGridMap, coarseMapBox, scanIndices,
                    boundingBox.mMin, x, y);

                /* Ignore the score of the low-resolution grid cell
                 * if the score is below a maximum score */
                if (normalizedScore <= scoreMax) {
                    /* Update the number of the ignored nodes */
                    numOfIgnoredNodes++;
                    continue;
                }

                /* Evaluate the score using the high-resolution grid map */
                /* Update the maximum score and search window index */
                this->EvaluateHighResolutionMap(
                    gridMap, boundingBox, scanIndices, Point2D<int>::Zero,
                    x, y, t, bestWinX, bestWinY, bestWinTheta, scoreMax);
                /* Update the number of the processed nodes */
                numOfProcessedNodes++;
            }
        }
    }

    /* The appropriate solution is found if the maximum score is
     * larger than (not larger than or equal to) the score threshold */
    const bool poseFound = scoreMax > normalizedScoreThreshold;
    /* Compute the best sensor pose */
    const RobotPose2D<double> bestSensorPose {
        sensorPose.mX + bestWinX * stepX,
        sensorPose.mY + bestWinY * stepY,
        sensorPose.mTheta + bestWinTheta * stepTheta };

    /* Compute the estimated robot pose and the likelihood value */
    estimatedPose = MoveBackward(
        bestSensorPose, scanData->RelativeSensorPose());
    likelihoodValue = this->mLikelihoodFunc->Likelihood(
        gridMap, scanData, bestSensorPose);
    normalizedLikelihood = likelihoodValue / scanData->NumOfScans();
    normalizedScore = scoreMax;

    /* Update the distribution and counter metrics */
    this->mMetrics.mDiffTranslation->Observe(
        Distance(initialPose, estimatedPose));
    this->mMetrics.mDiffRotation->Observe(
        std::abs(initialPose.mTheta - estimatedPose.mTheta));
    this->mMetrics.mWinSizeX->Observe(winX);
    this->mMetrics.mWinSizeY->Observe(winY);
    this->mMetrics.mWinSizeTheta->Observe(winTheta);
    this->mMetrics.mStepSizeX->Observe(stepX);
    this->mMetrics.mStepSizeY->Observe(stepY);
    this->mMetrics.mStepSizeTheta->Observe(stepTheta);
    this->mMetrics.mNumOfIgnoredNodes->Increment(numOfIgnoredNodes);
    this->mMetrics.mNumOfProcessedNodes->Increment(numOfProcessedNodes);
}

/* Compute the search step */
void ScanMatcherCorrelative::ComputeSearchStep(
    const GridMapInterface& gridMap,
    const Sensor::ScanDataPtr<double>& scanData,
    double& stepX,
    double& stepY,
    double& stepTheta) const
{
    /* Determine the search step */
    const auto maxRangeIt = std::max_element(
        scanData->Ranges().cbegin(), scanData->Ranges().cend());
    const double maxRange = std::min(*maxRangeIt, this->mScanRangeMax);
    const double theta = gridMap.Resolution() / maxRange;

    stepX = gridMap.Resolution();
    stepY = gridMap.Resolution();
    stepTheta = std::acos(1.0 - 0.5 * theta * theta);

    return;
}

/* Compute the grid cell indices for scan points */
void ScanMatcherCorrelative::ComputeScanIndices(
    const GridMapInterface& gridMap,
    const RobotPose2D<double>& sensorPose,
    const Sensor::ScanDataPtr<double>& scanData,
    std::vector<Point2D<int>>& scanIndices) const
{
    /* Compute the grid cell indices for scan points */
    scanIndices.clear();

    const std::size_t numOfScans = scanData->NumOfScans();

    for (std::size_t i = 0; i < numOfScans; ++i) {
        const double range = scanData->RangeAt(i);

        if (range >= this->mScanRangeMax)
            continue;

        Point2D<double> hitPoint = scanData->HitPoint(sensorPose, i);
        Point2D<int> hitIdx = gridMap.PositionToIndex(hitPoint.mX, hitPoint.mY);
        scanIndices.push_back(std::move(hitIdx));
    }

    return;
}

/* Compute the normalized scan matching score based on the
 * already projected scan points (indices) and index offsets */
double ScanMatcherCorrelative::ComputeScore(
    const GridMapInterface& gridMap,
    const BoundingBox<int>& boundingBox,
    const std::vector<Point2D<int>>& scanIndices,
    const Point2D<int>& scanIdxOffset,
    const int offsetX,
    const int offsetY) const
{
    /* Evaluate the matching score based on the occupancy probability value */
    const double unknownProb = gridMap.UnknownProbability();
    double sumScore = 0.0;

    for (const auto& hitIdx : scanIndices) {
        /* Compute the grid cell index based on the offsets */
        const Point2D<int> gridCellIdx {
            hitIdx.mX - scanIdxOffset.mX + offsetX,
            hitIdx.mY - scanIdxOffset.mY + offsetY };
        if (!boundingBox.IsInside(gridCellIdx.mX, gridCellIdx.mY))
            continue;

        /* Retrieve the occupancy probability value from the grid map */
        const double prob = gridMap.ProbabilityOr(
            gridCellIdx.mY, gridCellIdx.mX, unknownProb);
        if (prob == unknownProb)
            continue;

        /* Only grid cells that are observed at least once and that have known
         * occupancy probability values are considered in the computation */
        sumScore += prob;
    }

    /* Normalize the score function */
    const double normalizedScore = sumScore / scanIndices.size();
    return normalizedScore;
}

/* Evaluate the matching score using high-resolution grid map */
void ScanMatcherCorrelative::EvaluateHighResolutionMap(
    const GridMapInterface& gridMap,
    const BoundingBox<int>& boundingBox,
    const std::vector<Point2D<int>>& scanIndices,
    const Point2D<int>& scanIdxOffset,
    const int offsetX,
    const int offsetY,
    const int offsetTheta,
    int& maxWinX,
    int& maxWinY,
    int& maxWinTheta,
    double& maxScore) const
{
    /* Search inside the relatively small area */
    for (int x = offsetX; x < offsetX + this->mLowResolution; ++x) {
        for (int y = offsetY; y < offsetY + this->mLowResolution; ++y) {
            /* Evaluate the normalized matching score */
            const double normalizedScore = this->ComputeScore(
                gridMap, boundingBox, scanIndices, scanIdxOffset, x, y);

            /* Update the maximum score and search window index */
            if (maxScore < normalizedScore) {
                maxScore = normalizedScore;
                maxWinX = x;
                maxWinY = y;
                maxWinTheta = offsetTheta;
            }
        }
    }

    return;
}

} /* namespace Mapping */
} /* namespace MyGMapping */
