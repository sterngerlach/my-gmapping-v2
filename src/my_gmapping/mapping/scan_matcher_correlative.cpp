
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

    /* Register the value sequence metrics */
    this->mInputSetupTime = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".InputSetupTime");
    this->mOptimizationTime = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".OptimizationTime");

    this->mDiffTranslation = pMetricManager->AddValueSequence<float>(
        scanMatcherName + ".DiffTranslation");
    this->mDiffRotation = pMetricManager->AddValueSequence<float>(
        scanMatcherName + ".DiffRotation");
    this->mWinSizeX = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".WinSizeX");
    this->mWinSizeY = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".WinSizeY");
    this->mWinSizeTheta = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".WinSizeTheta");
    this->mStepSizeX = pMetricManager->AddValueSequence<float>(
        scanMatcherName + ".StepSizeX");
    this->mStepSizeY = pMetricManager->AddValueSequence<float>(
        scanMatcherName + ".StepSizeY");
    this->mStepSizeTheta = pMetricManager->AddValueSequence<float>(
        scanMatcherName + ".StepSizeTheta");

    this->mNumOfIgnoredNodes = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".NumOfIgnoredNodes");
    this->mNumOfProcessedNodes = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".NumOfProcessedNodes");

    this->mScoreValue = pMetricManager->AddValueSequence<float>(
        scanMatcherName + ".ScoreValue");
    this->mLikelihoodValue = pMetricManager->AddValueSequence<float>(
        scanMatcherName + ".LikelihoodValue");
    this->mNumOfScans = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".NumOfScans");
}

/* Reserve the buffer to store the processing times */
void ScanMatcherCorrelativeMetrics::Resize(
    const std::size_t numOfParticles)
{
    this->mTimes.resize(numOfParticles);
    this->mParams.resize(numOfParticles);
}

/* Set the processing times for each particle */
void ScanMatcherCorrelativeMetrics::SetTimes(
    const std::size_t idx, const Times& times)
{
    this->mTimes[idx] = times;
}

/* Set the parameter settings for each particle */
void ScanMatcherCorrelativeMetrics::SetParameters(
    const std::size_t idx, const Parameters& parameters)
{
    this->mParams[idx] = parameters;
}

/* Collect the particle-wise metrics and update the overall metrics */
void ScanMatcherCorrelativeMetrics::Update(
    const std::size_t bestIdx)
{
    auto collectTimes = [&](
        std::function<int(int, const Times&)> selector) {
        return std::accumulate(this->mTimes.begin(),
                               this->mTimes.end(), 0, selector); };
    auto collectParams = [&](
        std::function<int(int, const Parameters&)> selector) {
            return std::accumulate(this->mParams.begin(),
                                   this->mParams.end(), 0, selector); };

    /* Compute the sum of the processing times */
    this->mInputSetupTime->Observe(collectTimes(
        [](int value, const Times& times) {
            return value + times.mInputSetupTime; }));
    this->mOptimizationTime->Observe(collectTimes(
        [](int value, const Times& times) {
            return value + times.mOptimizationTime; }));

    /* Store the metric values of the best particle */
    const auto& bestParticle = this->mParams[bestIdx];
    this->mDiffTranslation->Observe(bestParticle.mDiffTranslation);
    this->mDiffRotation->Observe(bestParticle.mDiffRotation);
    this->mWinSizeX->Observe(bestParticle.mWinSizeX);
    this->mWinSizeY->Observe(bestParticle.mWinSizeY);
    this->mWinSizeTheta->Observe(bestParticle.mWinSizeTheta);
    this->mStepSizeX->Observe(bestParticle.mStepSizeX);
    this->mStepSizeY->Observe(bestParticle.mStepSizeY);
    this->mStepSizeTheta->Observe(bestParticle.mStepSizeTheta);
    this->mScoreValue->Observe(bestParticle.mScoreValue);
    this->mLikelihoodValue->Observe(bestParticle.mLikelihoodValue);
    this->mNumOfScans->Observe(bestParticle.mNumOfScans);

    /* Store the collected metric values */
    this->mNumOfIgnoredNodes->Observe(collectParams(
        [](int value, const Parameters& params) {
            return value + params.mNumOfIgnoredNodes; }));
    this->mNumOfProcessedNodes->Observe(collectParams(
        [](int value, const Parameters& params) {
            return value + params.mNumOfProcessedNodes; }));
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
ScanMatchingResultVector ScanMatcherCorrelative::OptimizePose(
    const ScanMatchingQueryVector& queries,
    const Sensor::ScanDataPtr<double>& scanData)
{
    ScanMatchingResultVector results;
    results.resize(queries.size());

    /* Resize the buffer to store the metrics for each particle */
    this->mMetrics.Resize(queries.size());

    /* Optimize the pose for each particle (parallelized using OpenMP) */
#pragma omp parallel for
    for (std::size_t i = 0; i < queries.size(); ++i) {
        /* Create the timer for each thread */
        Metric::Timer timer;
        /* Metric values for each particle */
        ScanMatcherMetrics::Times times;
        ScanMatcherMetrics::Parameters params;

        const auto& initialPose = queries[i].mInitialPose;
        const auto& gridMap = queries[i].mGridMap;
        const int mapCols = gridMap.Cols();
        const int mapRows = gridMap.Rows();

        /* Create the grid map to store the intermediate result and the final
         * result (low-resolution coarse grid map) */
        ConstMap intermediateMap {
            gridMap.Resolution(), gridMap.BlockSize(),
            this->mCroppedMapSizeY, this->mCroppedMapSizeX };
        ConstMap precompMap {
            gridMap.Resolution(), gridMap.BlockSize(),
            this->mCroppedMapSizeY, this->mCroppedMapSizeX };

        /* Determine the desired size of the grid map for scan matching */
        const int mapColsMax = this->mUseCroppedMap ?
            this->mCroppedMapSizeX : mapCols;
        const int mapRowsMax = this->mUseCroppedMap ?
            this->mCroppedMapSizeY : mapRows;

        /* Compute the center position of the cropped grid map */
        const Point2D<int> desiredCenterIdx = gridMap.PositionToIndex(
            initialPose.mX, initialPose.mY);
        /* Use std::clamp() here since the above `centerIdx` could be
         * out-of-bounds (the current particle is outside of the grid map) */
        const Point2D<int> possibleIdxMin {
            std::clamp(desiredCenterIdx.mX - mapColsMax, 0, mapCols - 1),
            std::clamp(desiredCenterIdx.mY - mapRowsMax, 0, mapRows - 1) };
        const Point2D<int> possibleIdxMax {
            std::clamp(desiredCenterIdx.mX + mapColsMax, 1, mapCols),
            std::clamp(desiredCenterIdx.mY + mapRowsMax, 1, mapRows) };
        const Point2D<int> centerIdx {
            (possibleIdxMax.mX + possibleIdxMin.mX) / 2,
            (possibleIdxMax.mY + possibleIdxMin.mY) / 2 };

        /* Crop the grid map around the center position */
        const Point2D<int> idxMin {
            std::max(centerIdx.mX - mapColsMax / 2, 0),
            std::max(centerIdx.mY - mapRowsMax / 2, 0) };
        const Point2D<int> idxMax {
            std::min(centerIdx.mX + mapColsMax / 2, mapCols),
            std::min(centerIdx.mY + mapRowsMax / 2, mapRows) };
        const BoundingBox<int> boundingBox { idxMin, idxMax };

        /* Precompute the coarser grid map from the cropped grid map */
        PrecomputeGridMap(gridMap, this->mLowResolution, idxMin,
                          intermediateMap, precompMap);

        /* Update the processing time for setting up the input */
        times.mInputSetupTime = timer.ElapsedMicro();
        timer.Start();

        /* Optimize the particle pose based on the correlative scan matching */
        results[i] = this->OptimizePoseCore(gridMap, boundingBox, precompMap,
                                            scanData, initialPose, 0.0, params);

        /* Update the processing time for the optimization */
        times.mOptimizationTime = timer.ElapsedMicro();
        timer.Stop();

        /* Set the metrics */
        this->mMetrics.SetTimes(i, times);
        this->mMetrics.SetParameters(i, params);
    }

    /* Determine the maximum likelihood value and its corresponding score */
    const auto bestIt = std::max_element(
        results.begin(), results.end(),
        [](const ScanMatchingResult& lhs, const ScanMatchingResult& rhs) {
            return lhs.mLikelihood < rhs.mLikelihood; });
    const auto bestIdx = std::distance(results.begin(), bestIt);

    /* Update the metrics */
    this->mMetrics.Update(bestIdx);

    return results;
}

/* Optimize the particle pose based on the correlative scan matching */
ScanMatchingResult ScanMatcherCorrelative::OptimizePoseCore(
    const GridMapInterface& gridMap,
    const BoundingBox<int>& boundingBox,
    const GridMapInterface& coarseGridMap,
    const Sensor::ScanDataPtr<double>& scanData,
    const RobotPose2D<double>& initialPose,
    const double normalizedScoreThreshold,
    ScanMatcherMetrics::Parameters& params)
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
    double scoreMax = normalizedScoreThreshold * scanData->NumOfScans();
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
                /* Evaluate the matching score
                 * Add index offset `boundingBox.mMin` since the grid cell
                 * indices for the scan points `scanIndices` are computed
                 * from the original grid map `gridMap` */
                const BoundingBox<int>& coarseMapBox {
                    0, 0, coarseGridMap.Cols(), coarseGridMap.Rows() };
                const double score = this->ComputeScore(
                    coarseGridMap, coarseMapBox, scanIndices,
                    boundingBox.mMin, x, y);

                /* Ignore the score of the low-resolution grid cell
                 * if the score is below a maximum score */
                if (score <= scoreMax) {
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

    /* Compute the best sensor pose */
    const RobotPose2D<double> bestSensorPose {
        sensorPose.mX + bestWinX * stepX,
        sensorPose.mY + bestWinY * stepY,
        sensorPose.mTheta + bestWinTheta * stepTheta };

    /* Compute the estimated robot pose and the likelihood value */
    const RobotPose2D<double> estimatedPose = MoveBackward(
        bestSensorPose, scanData->RelativeSensorPose());
    const RobotPose2D<double> diffPose =
        InverseCompound(initialPose, estimatedPose);

    /* Set the resulting likelihood */
    const double likelihood = this->mLikelihoodFunc->Likelihood(
        gridMap, scanData, bestSensorPose);
    const double normalizedLikelihood = likelihood / scanData->NumOfScans();
    /* Set the resulting score */
    const double score = scoreMax;
    const double normalizedScore = score / scanData->NumOfScans();

    /* Update the metrics */
    params.mDiffTranslation = Distance(diffPose);
    params.mDiffRotation = std::abs(diffPose.mTheta);
    params.mWinSizeX = winX;
    params.mWinSizeY = winY;
    params.mWinSizeTheta = winTheta;
    params.mStepSizeX = stepX;
    params.mStepSizeY = stepY;
    params.mStepSizeTheta = stepTheta;
    params.mNumOfIgnoredNodes = numOfIgnoredNodes;
    params.mNumOfProcessedNodes = numOfProcessedNodes;
    params.mScoreValue = normalizedScore;
    params.mLikelihoodValue = normalizedLikelihood;
    params.mNumOfScans = scanData->NumOfScans();

    return ScanMatchingResult { initialPose, estimatedPose,
                                normalizedLikelihood, likelihood,
                                normalizedScore, score };
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

/* Compute the scan matching score based on the already projected
 * scan points (indices) and index offsets */
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

    return sumScore;
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
            /* Evaluate the matching score */
            const double score = this->ComputeScore(
                gridMap, boundingBox, scanIndices, scanIdxOffset, x, y);

            /* Update the maximum score and search window index */
            if (maxScore < score) {
                maxScore = score;
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
