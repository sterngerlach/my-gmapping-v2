
/* grid_map_builder.cpp */

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "my_gmapping/util.hpp"
#include "my_gmapping/mapping/grid_map_builder.hpp"

namespace MyGMapping {
namespace Mapping {

/* Constructor */
GridMapBuilderMetrics::GridMapBuilderMetrics() :
    mInputScanDataCount(nullptr),
    mProcessCount(nullptr),
    mProcessTime(nullptr),
    mProcessScanTime(nullptr),
    mSamplingTime(nullptr),
    mScanDataSetupTime(nullptr),
    mScanMatchingTime(nullptr),
    mFinalScanMatchingTime(nullptr),
    mWeightUpdateTime(nullptr),
    mMapUpdateTime(nullptr),
    mLatestMapUpdateTime(nullptr),
    mResamplingTime(nullptr),
    mIntervalTravelDist(nullptr),
    mIntervalAngle(nullptr),
    mIntervalTime(nullptr),
    mNumOfScans(nullptr),
    mProcessFrame(nullptr),
    mEffectiveSampleSize(nullptr)
{
    /* Retrieve the metrics manager instance */
    auto* const pMetricManager = Metric::MetricManager::Instance();

    /* Register the counter metrics */
    this->mInputScanDataCount = pMetricManager->AddCounter(
        "GridMapBuilder.InputScanDataCount");
    this->mProcessCount = pMetricManager->AddCounter(
        "GridMapBuilder.ProcessCount");
    this->mProcessTime = pMetricManager->AddCounter(
        "GridMapBuilder.ProcessTime");

    /* Register the distribution metrics */
    this->mProcessScanTime = pMetricManager->AddDistribution(
        "GridMapBuilder.ProcessScanTime");
    this->mSamplingTime = pMetricManager->AddDistribution(
        "GridMapBuilder.SamplingTime");
    this->mScanDataSetupTime = pMetricManager->AddDistribution(
        "GridMapBuilder.ScanDataSetupTime");
    this->mScanMatchingTime = pMetricManager->AddDistribution(
        "GridMapBuilder.ScanMatchingTime");
    this->mFinalScanMatchingTime = pMetricManager->AddDistribution(
        "GridMapBuilder.FinalScanMatchingTime");
    this->mWeightUpdateTime = pMetricManager->AddDistribution(
        "GridMapBuilder.WeightUpdateTime");
    this->mMapUpdateTime = pMetricManager->AddDistribution(
        "GridMapBuilder.MapUpdateTime");
    this->mLatestMapUpdateTime = pMetricManager->AddDistribution(
        "GridMapBuilder.LatestMapUpdateTime");
    this->mResamplingTime = pMetricManager->AddDistribution(
        "GridMapBuilder.ResamplingTime");
    this->mIntervalTravelDist = pMetricManager->AddDistribution(
        "GridMapBuilder.IntervalTravelDist");
    this->mIntervalAngle = pMetricManager->AddDistribution(
        "GridMapBuilder.IntervalAngle");
    this->mIntervalTime = pMetricManager->AddDistribution(
        "GridMapBuilder.IntervalTime");
    this->mNumOfScans = pMetricManager->AddDistribution(
        "GridMapBuilder.NumOfScans");

    /* Register the value sequence metrics */
    this->mProcessFrame = pMetricManager->AddValueSequence<int>(
        "GridMapBuilder.ProcessFrame");
    this->mEffectiveSampleSize = pMetricManager->AddValueSequence<float>(
        "GridMapBuilder.EffectiveSampleSize");
}

/* Constructor */
GridMapBuilder::GridMapBuilder(
    std::unique_ptr<MotionModel>&& motionModel,
    std::unique_ptr<LikelihoodFunction>&& likelihoodFunc,
    std::unique_ptr<ScanMatcher>&& scanMatcher,
    std::unique_ptr<ScanMatcher>&& finalScanMatcher,
    std::unique_ptr<ScanInterpolator>&& scanInterpolator,
    std::unique_ptr<CovarianceEstimator>&& covarianceEstimator,
    MapBuilder&& mapBuilder,
    const WeightNormalizationType weightNormalizationType,
    const int numOfParticles,
    const bool useLatestMap,
    const std::size_t numOfScansForLatestMap,
    const RobotPose2D<double>& initialPose,
    const double mapCellSize,
    const double updateThresholdTravelDist,
    const double updateThresholdAngle,
    const double updateThresholdTime,
    const double resampleThreshold,
    const double degenerationThreshold) :
    mProcessCount(0),
    mMotionModel(std::move(motionModel)),
    mLikelihoodFunc(std::move(likelihoodFunc)),
    mScanMatcher(std::move(scanMatcher)),
    mFinalScanMatcher(std::move(finalScanMatcher)),
    mScanInterpolator(std::move(scanInterpolator)),
    mCovarianceEstimator(std::move(covarianceEstimator)),
    mMapBuilder(std::move(mapBuilder)),
    mWeightNormalizationType(weightNormalizationType),
    mRandEngine(std::random_device()()),
    mUseLatestMap(useLatestMap),
    mNumOfScansForLatestMap(numOfScansForLatestMap),
    mLastOdomPose(0.0, 0.0, 0.0),
    mAccumulatedTravelDist(0.0),
    mAccumulatedAngle(0.0),
    mLastMapUpdateOdomPose(0.0, 0.0, 0.0),
    mLastMapUpdateTime(0.0),
    mUpdateThresholdTravelDist(updateThresholdTravelDist),
    mUpdateThresholdAngle(updateThresholdAngle),
    mUpdateThresholdTime(updateThresholdTime),
    mResampleThreshold(resampleThreshold),
    mDegenerationThreshold(degenerationThreshold)
{
    XAssert(numOfParticles > 0,
            "The number of particles must be positive");
    XAssert(numOfScansForLatestMap > 0,
            "The number of scans used for the latest map must be positive");
    XAssert(mapCellSize > 0.0,
            "The resolution of the grid map must be positive");
    XAssert(updateThresholdTravelDist >= 0.0,
            "The accumulated travel distance used as a threshold for "
            "executing the scan matching and map update must be "
            "greater than or equal to 0");
    XAssert(updateThresholdAngle >= 0.0,
            "The accumulated angle used as a threshold for "
            "executing the scan matching and map update must be "
            "greater than or equal to 0");
    XAssert(updateThresholdTime >= 0.0,
            "The elapsed time since the last map update used as a threshold for "
            "executing the scan matching and map update must be "
            "greater than or equal to 0");
    XAssert(resampleThreshold >= 0.0 && resampleThreshold <= 1.0,
            "The threshold for particle resampling must be within the range "
            "between 0.0 and 1.0");

    /* Initialize particles */
    this->mParticles.clear();
    this->mParticles.reserve(numOfParticles);

    /* Create root trajectory node and initial map */
    GridMapType initialMap { mapCellSize, 0.0, 0.0, 0.0, 0.0 };
    GridMapType initialLatestMap { mapCellSize, 0.0, 0.0, 0.0, 0.0 };
    const double initialWeight = 1.0 / static_cast<double>(numOfParticles);

    /* Insert particles */
    for (int i = 0; i < numOfParticles; ++i)
        this->mParticles.emplace_back(initialMap, initialLatestMap,
                                      initialWeight, initialPose, nullptr);
}

/* Process scan data */
bool GridMapBuilder::ProcessScan(
    const Sensor::ScanDataPtr<double>& rawScanData,
    const RobotPose2D<double>& odomPose)
{
    /* Create the timer to collect metrics */
    Metric::Timer outerTimer;
    Metric::Timer timer;

    /* Calculate the relative odometry pose */
    const RobotPose2D<double> relOdomPose = (this->mProcessCount == 0) ?
        RobotPose2D<double>(0.0, 0.0, 0.0) :
        InverseCompound(this->mLastOdomPose, odomPose);

    /* Set the previous odometry pose */
    this->mLastOdomPose = odomPose;

    /* Accumulate the linear and angular travel distance */
    this->mAccumulatedTravelDist += Distance(relOdomPose);
    this->mAccumulatedAngle += std::fabs(relOdomPose.mTheta);

    /* Compute the elapsed time since the last map update */
    const double elapsedTime = (this->mProcessCount == 0) ?
        0.0 : rawScanData->TimeStamp() - this->mLastMapUpdateTime;

    /* Particle map is updated only when the particle moved more than
     * the specified distance or the specified time has elapsed
     * since the last map update */
    const bool travelDistThreshold =
        this->mAccumulatedTravelDist >= this->mUpdateThresholdTravelDist;
    const bool angleThreshold =
        this->mAccumulatedAngle >= this->mUpdateThresholdAngle;
    const bool timeThreshold =
        elapsedTime >= this->mUpdateThresholdTime;
    const bool isFirstScan = this->mProcessCount == 0;
    const bool mapUpdateNeeded =
        (travelDistThreshold || angleThreshold ||
         timeThreshold || isFirstScan) && (elapsedTime >= 0.0);

    /* Update the total number of the input scan data */
    this->mMetrics.mInputScanDataCount->Increment();

    if (!mapUpdateNeeded) {
        /* Update the metrics */
        this->mMetrics.mProcessTime->Increment(outerTimer.ElapsedMicro());
        return false;
    }

    /* Compute the relative odometry pose since the last map update */
    const RobotPose2D<double> relPoseFromLastMapUpdate =
        InverseCompound(this->mLastMapUpdateOdomPose, odomPose);

    /* Restart the timer */
    timer.Start();

    /* Sample the particle poses from the odometry data
     * Update the current particle pose */
    for (auto& currentParticle : this->mParticles)
        currentParticle.Pose() = this->mMotionModel->SamplePose(
            currentParticle.Pose(), relPoseFromLastMapUpdate);

    /* Update the minimum position of the grid map */
    if (isFirstScan)
        for (auto& currentParticle : this->mParticles)
            currentParticle.Map().SetMinPos(
                currentParticle.Pose().mX, currentParticle.Pose().mY);

    /* Update the total processing time for the particle sampling */
    this->mMetrics.mSamplingTime->Observe(timer.ElapsedMicro());
    timer.Start();

    /* Interpolate scan data if necessary */
    auto scanData = (this->mScanInterpolator != nullptr) ?
        this->mScanInterpolator->Interpolate(rawScanData) : rawScanData;

    /* Update the total processing time for setting up the scan data */
    this->mMetrics.mScanDataSetupTime->Observe(timer.ElapsedMicro());
    timer.Start();

    /* Update the particle map with the latest scan data */
    if (this->mUseLatestMap)
        this->UpdateLatestMaps(this->mLatestScanData);

    /* Update the total processing time for updating the latest map */
    this->mMetrics.mLatestMapUpdateTime->Observe(timer.ElapsedMicro());
    timer.Stop();

    /* Improve particle poses by scan matching */
    this->ExecuteScanMatching(scanData);

    /* Start the timer */
    timer.Start();

    /* Update the particle map */
    this->UpdateGridMaps(scanData);
    /* Update the total processing time for updating the grid map */
    this->mMetrics.mMapUpdateTime->Observe(timer.ElapsedMicro());
    timer.Start();

    /* Resample if necessary */
    this->ResampleParticles();
    /* Update the total processing time for the particle resampling */
    this->mMetrics.mResamplingTime->Observe(timer.ElapsedMicro());
    timer.Stop();

    /* Update the deque that stores the latest scan data */
    if (this->mUseLatestMap)
        this->UpdateLatestScans(scanData);

    /* Update the metrics */
    this->mMetrics.mProcessCount->Increment();
    this->mMetrics.mProcessTime->Increment(outerTimer.ElapsedMicro());
    this->mMetrics.mProcessScanTime->Observe(outerTimer.ElapsedMicro());
    this->mMetrics.mIntervalTravelDist->Observe(this->mAccumulatedTravelDist);
    this->mMetrics.mIntervalAngle->Observe(this->mAccumulatedAngle);
    this->mMetrics.mIntervalTime->Observe(elapsedTime);
    this->mMetrics.mNumOfScans->Observe(scanData->NumOfScans());
    this->mMetrics.mProcessFrame->Observe(
        this->mMetrics.mInputScanDataCount->Value() - 1);

    /* Update miscellaneous parameters */
    this->mProcessCount += 1;
    this->mAccumulatedTravelDist = 0.0;
    this->mAccumulatedAngle = 0.0;
    this->mLastMapUpdateOdomPose = odomPose;
    this->mLastMapUpdateTime = scanData->TimeStamp();

    std::cerr << "Processing frame: " << this->mProcessCount << std::endl;

    return true;
}

/* Get the index of the best particle */
std::size_t GridMapBuilder::BestParticleIndex() const
{
    const auto bestIter = std::max_element(
        this->mParticles.cbegin(),
        this->mParticles.cend(),
        [](const Particle& lhs, const Particle& rhs) {
            return lhs.Weight() < rhs.Weight(); });
    const std::size_t bestIdx = static_cast<std::size_t>(
        std::distance(std::cbegin(this->mParticles), bestIter));

    return bestIdx;
}

/* Get the estimated trajectory of the specified particle */
std::vector<TimeStampedPose> GridMapBuilder::ParticleTrajectoryWithTimeStamp(
    std::size_t particleIdx) const
{
    const Particle& particle = this->mParticles[particleIdx];
    std::shared_ptr<TrajectoryNode> nodePtr = particle.Node();

    std::vector<TimeStampedPose> particleTrajectory;

    while (nodePtr != nullptr) {
        particleTrajectory.push_back(nodePtr->StampedPose());
        nodePtr = nodePtr->Parent();
    }

    /* Reverse the trajectory so that the first particle pose comes at first */
    std::reverse(std::begin(particleTrajectory),
                 std::end(particleTrajectory));

    return particleTrajectory;
}

/* Get the estimated trajectory of the specified particle */
std::vector<RobotPose2D<double>> GridMapBuilder::ParticleTrajectory(
    std::size_t particleIdx) const
{
    const Particle& particle = this->mParticles[particleIdx];
    std::shared_ptr<TrajectoryNode> nodePtr = particle.Node();

    std::vector<RobotPose2D<double>> particleTrajectory;

    while (nodePtr != nullptr) {
        particleTrajectory.emplace_back(nodePtr->Pose());
        nodePtr = nodePtr->Parent();
    }

    /* Reverse the trajectory so that the first particle pose comes at first */
    std::reverse(std::begin(particleTrajectory),
                 std::end(particleTrajectory));

    return particleTrajectory;
}

/* Update the deque that stores the latest scan data */
void GridMapBuilder::UpdateLatestScans(
    const Sensor::ScanDataPtr<double>& scanData)
{
    /* Push the latest scan data to the head of the deque
     * tail of the deque is the oldest */
    this->mLatestScanData.push_front(scanData);

    if (this->mLatestScanData.size() > this->mNumOfScansForLatestMap)
        this->mLatestScanData.pop_back();
}

/* Execute scan matching and update particle weights */
void GridMapBuilder::ExecuteScanMatching(
    const Sensor::ScanDataPtr<double>& scanData)
{
    /* Create the timer */
    Metric::Timer timer;

    const std::size_t numOfParticles = this->mParticles.size();

    /* Setup scan matching outputs */
    std::vector<double> particleWeights;
    particleWeights.reserve(numOfParticles);
    std::vector<RobotPose2D<double>> estimatedPoses;
    estimatedPoses.reserve(numOfParticles);

    /* Setup scan matching inputs */
    std::vector<const GridMapType*> particleMaps;
    particleMaps.reserve(numOfParticles);
    std::vector<RobotPose2D<double>> initialPoses;
    initialPoses.reserve(numOfParticles);

    for (std::size_t i = 0; i < numOfParticles; ++i) {
        /* Use the grid map constructed from the latest scans
         * instead of the entire grid map */
        if (this->mUseLatestMap)
            particleMaps.push_back(&this->mParticles[i].LatestMap());
        else
            particleMaps.push_back(&this->mParticles[i].Map());

        initialPoses.push_back(this->mParticles[i].Pose());
    }

    /* Execute scan matching for particles */
    this->mScanMatcher->OptimizePose(
        numOfParticles, particleMaps, scanData,
        initialPoses, estimatedPoses, particleWeights);

    /* Update the total processing time for the scan matching */
    this->mMetrics.mScanMatchingTime->Observe(timer.ElapsedMicro());
    timer.Start();

    /* Update the initial poses */
    std::vector<RobotPose2D<double>> refinedInitialPoses;
    refinedInitialPoses.reserve(numOfParticles);

    for (std::size_t i = 0; i < numOfParticles; ++i)
        refinedInitialPoses.push_back(estimatedPoses[i]);

    /* Execute the final scan matching for particles */
    this->mFinalScanMatcher->OptimizePose(
        numOfParticles, particleMaps, scanData,
        refinedInitialPoses, estimatedPoses, particleWeights);

    /* Update the particle pose and trajectory node */
    for (std::size_t i = 0; i < numOfParticles; ++i) {
        /* Check the degeneration using the pose covariance */
        const Eigen::Matrix3d covarianceMat =
            this->mCovarianceEstimator->ComputeCovariance(
                this->mParticles[i].Map(), scanData, estimatedPoses[i]);
        const bool degenerationDetected =
            this->CheckDegeneration(covarianceMat);

        /* Update the particle pose if there is no degeneration */
        if (!degenerationDetected)
            this->mParticles[i].Pose() = estimatedPoses[i];

        /* Update the trajectory node */
        auto pNewNode = std::make_shared<TrajectoryNode>(
            this->mParticles[i].Node(), this->mParticles[i].Pose(),
            scanData->TimeStamp());
        this->mParticles[i].Node() = pNewNode;
    }

    /* Update the total processing time for the final scan matching */
    this->mMetrics.mFinalScanMatchingTime->Observe(timer.ElapsedMicro());
    timer.Start();

    /* Compute the log likelihood using the entire grid map if necessary */
    if (this->mUseLatestMap) {
        for (std::size_t i = 0; i < numOfParticles; ++i) {
            /* Compute the sensor pose from the particle pose */
            const RobotPose2D<double> sensorPose =
                Compound(estimatedPoses[i], scanData->RelativeSensorPose());

            /* Compute the log likelihood */
            const double logLikelihood = this->mLikelihoodFunc->Likelihood(
                this->mParticles[i].Map(), scanData, sensorPose);

            /* Update the log likelihood */
            particleWeights[i] = logLikelihood;
        }
    }

    /* Normalize the particle weights */
    this->mWeightNormalizer.NormalizeWeights(
        this->mWeightNormalizationType, particleWeights);

    /* Update the particle weights */
    for (std::size_t i = 0; i < this->mParticles.size(); ++i)
        this->mParticles[i].SetWeight(particleWeights[i]);

    /* Update the total processing time for updating the particle weights */
    this->mMetrics.mWeightUpdateTime->Observe(timer.ElapsedMicro());
    timer.Stop();
}

/* Integrate the scan data to the particle map */
void GridMapBuilder::UpdateGridMap(
    Particle& particle,
    const Sensor::ScanDataPtr<double>& scanData)
{
    /* Just forward to the map builder */
    this->mMapBuilder.UpdateGridMap(
        particle.Map(), particle.Pose(), scanData);
}

/* Update the grid maps for all particles */
void GridMapBuilder::UpdateGridMaps(
    const Sensor::ScanDataPtr<double>& scanData)
{
    /* Update the particle map */
    for (auto& currentParticle : this->mParticles)
        this->UpdateGridMap(currentParticle, scanData);
}

/* Update the particle map with the multiple latest scans */
void GridMapBuilder::UpdateLatestMap(
    Particle& particle,
    const ScanDataDeque& latestScanData)
{
    /* Return if the scan data is unavailable (first iteration) */
    if (latestScanData.empty())
        return;

    auto pNode = particle.Node();
    std::deque<RobotPose2D<double>> latestPoses;

    for (std::size_t i = 0; i < latestScanData.size(); ++i) {
        latestPoses.push_front(pNode->Pose());
        pNode = pNode->Parent();
    }

    /* Just forward to the map builder */
    this->mMapBuilder.UpdateLatestMap(
        particle.LatestMap(), latestPoses, latestScanData);
}

/* Update the latest maps for all particles */
void GridMapBuilder::UpdateLatestMaps(
    const ScanDataDeque& latestScanData)
{
    for (auto& currentParticle : this->mParticles)
        this->UpdateLatestMap(currentParticle, latestScanData);
}

/* Resample particles according to their weights if necessary */
void GridMapBuilder::ResampleParticles()
{
    static std::uniform_real_distribution uniformDist { 0.0, 1.0 };

    const std::size_t numOfParticles = this->mParticles.size();
    const double resampleThreshold =
        static_cast<double>(this->mParticles.size()) *
        this->mResampleThreshold;
    const double numOfEffectiveSamples =
        this->CalculateEffectiveSampleSize();

    /* Update the effective sample size */
    this->mMetrics.mEffectiveSampleSize->Observe(numOfEffectiveSamples);

    /* Resample only when the effective sample size becomes
     * less than the threshold */
    if (numOfEffectiveSamples >= resampleThreshold)
        return;

    std::cerr << std::fixed << std::setprecision(2);
    std::cerr << "Resample (Effective Sample Size: "
              << numOfEffectiveSamples << ")" << std::endl;

    /* Calculate the width for low-variance sampling */
    const double weightSum = std::accumulate(
        this->mParticles.begin(), this->mParticles.end(), 0.0,
        [](double accWeight, const Particle& particle) {
            return accWeight + particle.Weight(); });
    const double weightUnit =
        weightSum / static_cast<double>(numOfParticles);

    /* Perform low-variance sampling */
    std::vector<int> resampledIndices;
    resampledIndices.reserve(numOfParticles);

    double accWeight = uniformDist(this->mRandEngine) * weightUnit;
    double accParticleWeight = this->mParticles[0].Weight();
    int particleIdx = 0;

    for (int i = 0; i < static_cast<int>(numOfParticles); ++i) {
        while (accWeight > accParticleWeight)
            accParticleWeight += this->mParticles[++particleIdx].Weight();

        resampledIndices.push_back(particleIdx);
        accWeight += weightUnit;
    }

    /* Replace particles */
    for (int i = static_cast<int>(numOfParticles - 1); i >= 0; --i)
        if (i != resampledIndices[i])
            this->mParticles[i] = this->mParticles[resampledIndices[i]];
}

/* Calculate the effective sample size of the particles */
double GridMapBuilder::CalculateEffectiveSampleSize() const
{
    const double sumOfSquaredWeights = std::accumulate(
        this->mParticles.cbegin(), this->mParticles.cend(), 0.0,
        [](double accWeight, const Particle& particle) {
            return accWeight + particle.Weight() * particle.Weight(); });
    return 1.0 / sumOfSquaredWeights;
}

/* Check the degeneration */
bool GridMapBuilder::CheckDegeneration(
    const Eigen::Matrix3d& poseCovarianceMat) const
{
    /* Compute the eigenvalues using the pose covariance matrix */
    const Eigen::Vector2cd covEigenvalues =
        poseCovarianceMat.block<2, 2>(0, 0).eigenvalues();
    const double minEigen = std::min(
        covEigenvalues[0].real(), covEigenvalues[1].real());
    const double maxEigen = std::max(
        covEigenvalues[0].real(), covEigenvalues[1].real());
    /* Compute the ratio between two eigenvalues */
    const double eigenRatio = maxEigen / minEigen;
    /* Degeneration is detected if the ratio is considerably large */
    return eigenRatio > this->mDegenerationThreshold;
}

/*
 * Utility function implementations
 */

/* Compute the maximum of a 'winSize' pixel wide row at each pixel */
void SlidingWindowMaxRow(const GridMapType& gridMap,
                         const int winSize,
                         const Point2D<int>& gridCellIdxMin,
                         ConstMapType& intermediateMap)
{
    /* Each cell in the grid map stores an occupancy probability value
     * using the data type `StorageType` (defined as std::uint16_t) */
    using StorageType = GridMapType::StorageType;

    /* Compute the maximum for each column */
    const StorageType unknownRawValue = gridMap.UnknownRawValue();
    int colIdx = 0;

    std::function<StorageType(int)> inFunc =
        [&colIdx, &gridMap, &gridCellIdxMin, unknownRawValue](int rowIdx) {
        return gridMap.RawValue(gridCellIdxMin.mX + colIdx,
                                gridCellIdxMin.mY + rowIdx,
                                unknownRawValue);
    };

    std::function<void(int, StorageType)> outFunc =
        [&colIdx, &gridMap, &gridCellIdxMin, &intermediateMap, unknownRawValue]
        (int rowIdx, StorageType maxValue) {
        const Point2D<int> patchIdx = gridMap.CellIndexToPatchIndex(
            gridCellIdxMin.mX + colIdx, gridCellIdxMin.mY + rowIdx);
        const bool isAllocated = gridMap.PatchIsAllocated(patchIdx);

        if (isAllocated)
            intermediateMap.SetRawValue(colIdx, rowIdx, maxValue);
    };

    const int numOfGridCellsX = intermediateMap.NumCellsX();
    const int numOfGridCellsY = intermediateMap.NumCellsY();

    /* Apply the sliding window maximum function */
    for (colIdx = 0; colIdx < numOfGridCellsX; ++colIdx)
        SlidingWindowMax(inFunc, outFunc, numOfGridCellsY, winSize);
}

/* Compute the maximum of a 'winSize' pixel wide column at each pixel */
void SlidingWindowMaxCol(const ConstMapType& intermediateMap,
                         const int winSize,
                         ConstMapType& precompMap)
{
    /* Each cell in the grid map stores an occupancy probability value
     * using the data type `StorageType` (defined as std::uint16_t) */
    using StorageType = GridMapType::StorageType;

    /* Compute the maximum for each row */
    const StorageType unknownRawValue = intermediateMap.UnknownRawValue();
    int rowIdx = 0;

    std::function<StorageType(int)> inFunc =
        [&rowIdx, &intermediateMap, unknownRawValue](int colIdx) {
        return intermediateMap.RawValue(colIdx, rowIdx, unknownRawValue);
    };

    std::function<void(int, StorageType)> outFunc =
        [&rowIdx, &intermediateMap, &precompMap, unknownRawValue]
        (int colIdx, StorageType maxValue) {
        const Point2D<int> patchIdx =
            intermediateMap.CellIndexToPatchIndex(colIdx, rowIdx);
        const bool isAllocated = intermediateMap.PatchIsAllocated(patchIdx);

        if (isAllocated)
            precompMap.SetRawValue(colIdx, rowIdx, maxValue);
    };

    const int numOfGridCellsX = intermediateMap.NumCellsX();
    const int numOfGridCellsY = intermediateMap.NumCellsY();

    /* Apply the sliding window maximum function */
    for (rowIdx = 0; rowIdx < numOfGridCellsY; ++rowIdx)
        SlidingWindowMax(inFunc, outFunc, numOfGridCellsX, winSize);
}

/* Precompute grid map for efficiency */
void PrecomputeGridMap(const GridMapType& gridMap,
                       const int winSize,
                       const Point2D<int>& gridCellIdxMin,
                       ConstMapType& intermediateMap,
                       ConstMapType& precompMap)
{
    /* Check if the resulting grid map size is same as the intermediate map */
    Assert(intermediateMap.NumCellsX() == precompMap.NumCellsX());
    Assert(intermediateMap.NumCellsY() == precompMap.NumCellsY());
    /* Check if the minimum grid cell index is inside the grid map */
    Assert(gridCellIdxMin.mX < gridMap.NumCellsX());
    Assert(gridCellIdxMin.mY < gridMap.NumCellsY());
    /* Check if the maximum grid cell index is inside the grid map */
    Assert(gridCellIdxMin.mX + precompMap.NumCellsX() > 0);
    Assert(gridCellIdxMin.mY + precompMap.NumCellsY() > 0);

    /* Each pixel in the `precompMap` stores the maximum of the
     * occupancy probability values of 'winSize' * 'winSize'
     * box of pixels beginning there */

    /* Store the maximum of the 'winSize' pixel wide row */
    SlidingWindowMaxRow(gridMap, winSize, gridCellIdxMin, intermediateMap);
    /* Store the maximum of the 'winSize' pixel wide column */
    SlidingWindowMaxCol(intermediateMap, winSize, precompMap);

    /* Set the minimum position of the `precompMap` and `intermediateMap`
     * Note that the GridMap<T>::CellIndexToMapCoordinate() returns the
     * minimum position of the given grid cell */
    const Point2D<double> minPos =
        gridMap.CellIndexToMapCoordinate(gridCellIdxMin);
    intermediateMap.SetMinPos(minPos);
    precompMap.SetMinPos(minPos);
}

} /* namespace Mapping */
} /* namespace MyGMapping */
