
/* grid_map_builder.cpp */

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <unordered_set>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "my_gmapping/memory_usage.hpp"
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
    mEffectiveSampleSize(nullptr),
    mPhysicalMemoryUsage(nullptr),
    mGridMapMemoryUsage(nullptr),
    mLatestMapMemoryUsage(nullptr),
    mTrajectoryMemoryUsage(nullptr)
{
    /* Retrieve the metrics manager instance */
    auto* const pMetricManager = Metric::MetricManager::Instance();

    /* Register the counter metrics */
    this->mInputScanDataCount = pMetricManager->AddCounter(
        "GridMapBuilder.InputScanDataCount");
    this->mProcessCount = pMetricManager->AddCounter(
        "GridMapBuilder.ProcessCount");

    /* Register the value sequence metrics */
    this->mProcessTime = pMetricManager->AddValueSequence<int>(
        "GridMapBuilder.ProcessTime");
    this->mProcessScanTime = pMetricManager->AddValueSequence<int>(
        "GridMapBuilder.ProcessScanTime");
    this->mSamplingTime = pMetricManager->AddValueSequence<int>(
        "GridMapBuilder.SamplingTime");
    this->mScanDataSetupTime = pMetricManager->AddValueSequence<int>(
        "GridMapBuilder.ScanDataSetupTime");
    this->mScanMatchingTime = pMetricManager->AddValueSequence<int>(
        "GridMapBuilder.ScanMatchingTime");
    this->mFinalScanMatchingTime = pMetricManager->AddValueSequence<int>(
        "GridMapBuilder.FinalScanMatchingTime");
    this->mWeightUpdateTime = pMetricManager->AddValueSequence<int>(
        "GridMapBuilder.WeightUpdateTime");
    this->mMapUpdateTime = pMetricManager->AddValueSequence<int>(
        "GridMapBuilder.MapUpdateTime");
    this->mLatestMapUpdateTime = pMetricManager->AddValueSequence<int>(
        "GridMapBuilder.LatestMapUpdateTime");
    this->mResamplingTime = pMetricManager->AddValueSequence<int>(
        "GridMapBuilder.ResamplingTime");
    this->mIntervalTravelDist = pMetricManager->AddValueSequence<float>(
        "GridMapBuilder.IntervalTravelDist");
    this->mIntervalAngle = pMetricManager->AddValueSequence<float>(
        "GridMapBuilder.IntervalAngle");
    this->mIntervalTime = pMetricManager->AddValueSequence<float>(
        "GridMapBuilder.IntervalTime");
    this->mNumOfScans = pMetricManager->AddValueSequence<int>(
        "GridMapBuilder.NumOfScans");

    this->mProcessFrame = pMetricManager->AddValueSequence<int>(
        "GridMapBuilder.ProcessFrame");
    this->mEffectiveSampleSize = pMetricManager->AddValueSequence<float>(
        "GridMapBuilder.EffectiveSampleSize");

    this->mPhysicalMemoryUsage =
        pMetricManager->AddValueSequence<std::uint64_t>(
            "GridMapBuilder.PhysicalMemoryUsage");
    this->mGridMapMemoryUsage =
        pMetricManager->AddValueSequence<std::uint64_t>(
            "GridMapBuilder.GridMapMemoryUsage");
    this->mLatestMapMemoryUsage =
        pMetricManager->AddValueSequence<std::uint64_t>(
            "GridMapBuilder.LatestMapMemoryUsage");
    this->mTrajectoryMemoryUsage =
        pMetricManager->AddValueSequence<std::uint64_t>(
            "GridMapBuilder.TrajectoryMemoryUsage");
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
    mResolution(mapCellSize),
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
    GridMap initialMap;
    GridMap initialLatestMap;
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
        this->mMetrics.mProcessTime->Observe(outerTimer.ElapsedMicro());
        return false;
    }

    /* Compute the relative odometry pose since the last map update */
    const RobotPose2D<double> relPoseFromLastMapUpdate =
        InverseCompound(this->mLastMapUpdateOdomPose, odomPose);

    /* Restart the timer */
    timer.Start();

    /* Sample the particle poses from the odometry data
     * Update the current particle pose */
    for (auto& particle : this->mParticles)
        particle.mPose = this->mMotionModel->SamplePose(
            particle.mPose, relPoseFromLastMapUpdate);

    /* Initialize the grid map and adjust the positional offset */
    if (isFirstScan) {
        for (auto& particle : this->mParticles) {
            const int blockSize = 1 << 4;
            particle.mMap.Initialize(
                this->mResolution, blockSize, blockSize, blockSize,
                particle.mPose.mX, particle.mPose.mY);
            particle.mLatestMap.Initialize(
                this->mResolution, blockSize, blockSize, blockSize,
                particle.mPose.mX, particle.mPose.mY);
        }
    }

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
    this->mMetrics.mProcessTime->Observe(outerTimer.ElapsedMicro());
    this->mMetrics.mProcessScanTime->Observe(outerTimer.ElapsedMicro());
    this->mMetrics.mIntervalTravelDist->Observe(this->mAccumulatedTravelDist);
    this->mMetrics.mIntervalAngle->Observe(this->mAccumulatedAngle);
    this->mMetrics.mIntervalTime->Observe(elapsedTime);
    this->mMetrics.mNumOfScans->Observe(scanData->NumOfScans());
    this->mMetrics.mProcessFrame->Observe(
        this->mMetrics.mInputScanDataCount->Value() - 1);
    this->mMetrics.mPhysicalMemoryUsage->Observe(GetPhysicalMemoryUsage());

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
            return lhs.mWeight < rhs.mWeight; });
    const std::size_t bestIdx = static_cast<std::size_t>(
        std::distance(std::cbegin(this->mParticles), bestIter));

    return bestIdx;
}

/* Get the estimated trajectory of the specified particle */
std::vector<TimeStampedPose> GridMapBuilder::ParticleTrajectoryWithTimeStamp(
    std::size_t particleIdx) const
{
    const Particle& particle = this->mParticles[particleIdx];
    std::shared_ptr<TrajectoryNode> nodePtr = particle.mNode;

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
    std::shared_ptr<TrajectoryNode> nodePtr = particle.mNode;

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

    /* Setup the scan matching queries */
    ScanMatchingQueryVector queries;

    /* Use the grid map constructed from the latest scans
     * instead of the entire grid map if specified */
    for (std::size_t i = 0; i < numOfParticles; ++i)
        if (this->mUseLatestMap)
            queries.emplace_back(this->mParticles[i].mLatestMap,
                                 this->mParticles[i].mPose);
        else
            queries.emplace_back(this->mParticles[i].mMap,
                                 this->mParticles[i].mPose);

    /* Execute scan matching for particles */
    const ScanMatchingResultVector results =
        this->mScanMatcher->OptimizePose(queries, scanData);

    /* Update the total processing time for the scan matching */
    this->mMetrics.mScanMatchingTime->Observe(timer.ElapsedMicro());
    timer.Start();

    /* Setup the scan matching queries */
    ScanMatchingQueryVector finalQueries;

    for (std::size_t i = 0; i < numOfParticles; ++i)
        if (this->mUseLatestMap)
            finalQueries.emplace_back(this->mParticles[i].mLatestMap,
                                      results[i].mEstimatedPose);
        else
            finalQueries.emplace_back(this->mParticles[i].mMap,
                                      results[i].mEstimatedPose);

    /* Execute the final scan matching for particles */
    const ScanMatchingResultVector finalResults =
        this->mFinalScanMatcher->OptimizePose(finalQueries, scanData);

    /* Update the particle pose and trajectory node */
    for (std::size_t i = 0; i < numOfParticles; ++i) {
        /* Check the degeneration using the pose covariance */
        const RobotPose2D<double> sensorPose = Compound(
            finalResults[i].mEstimatedPose, scanData->RelativeSensorPose());
        const Eigen::Matrix3d covarianceMat =
            this->mCovarianceEstimator->ComputeCovariance(
                this->mParticles[i].mMap, scanData, sensorPose);
        const bool degenerationDetected =
            this->CheckDegeneration(covarianceMat);

        /* Update the particle pose if there is no degeneration */
        if (!degenerationDetected)
            this->mParticles[i].mPose = finalResults[i].mEstimatedPose;

        /* Update the trajectory node */
        auto pNewNode = std::make_shared<TrajectoryNode>(
            this->mParticles[i].mNode, this->mParticles[i].mPose,
            scanData->TimeStamp());
        this->mParticles[i].mNode = pNewNode;
    }

    /* Update the total processing time for the final scan matching */
    this->mMetrics.mFinalScanMatchingTime->Observe(timer.ElapsedMicro());
    timer.Start();

    /* Compute the log likelihood using the entire grid map if necessary */
    std::vector<double> particleWeights;
    particleWeights.resize(numOfParticles);

    for (std::size_t i = 0; i < numOfParticles; ++i) {
        if (this->mUseLatestMap) {
            /* Compute the sensor pose from the particle pose */
            const RobotPose2D<double> sensorPose = Compound(
                finalResults[i].mEstimatedPose,
                scanData->RelativeSensorPose());
            /* Compute the log likelihood */
            particleWeights[i] = this->mLikelihoodFunc->Likelihood(
                this->mParticles[i].mMap, scanData, sensorPose);
        } else {
            particleWeights[i] = finalResults[i].mLikelihood;
        }
    }

    /* Normalize the particle weights */
    this->mWeightNormalizer.NormalizeWeights(
        this->mWeightNormalizationType, particleWeights);

    /* Update the particle weights */
    for (std::size_t i = 0; i < numOfParticles; ++i)
        this->mParticles[i].mWeight = particleWeights[i];

    /* Update the total processing time for updating the particle weights */
    this->mMetrics.mWeightUpdateTime->Observe(timer.ElapsedMicro());
    timer.Stop();

    /* Update the memory usage for the trajectory nodes */
    this->mMetrics.mTrajectoryMemoryUsage->Observe(
        this->InspectTrajectoryMemoryUsage());
}

/* Update the grid maps for all particles */
void GridMapBuilder::UpdateGridMaps(
    const Sensor::ScanDataPtr<double>& scanData)
{
    /* Update the particle map */
    for (auto& particle : this->mParticles)
        this->mMapBuilder.UpdateGridMap(
            particle.mMap, particle.mPose, scanData);

    /* Inspect the memory usage for the grid maps */
    this->mMetrics.mGridMapMemoryUsage->Observe(
        this->InspectGridMapMemoryUsage());
}

/* Update the latest maps for all particles */
void GridMapBuilder::UpdateLatestMaps(
    const ScanDataDeque& latestScanData)
{
    if (latestScanData.empty())
        return;

    for (auto& particle : this->mParticles) {
        auto pNode = particle.mNode;
        std::deque<RobotPose2D<double>> latestPoses;

        for (std::size_t i = 0; i < latestScanData.size(); ++i) {
            latestPoses.push_front(pNode->Pose());
            pNode = pNode->Parent();
        }

        this->mMapBuilder.UpdateLatestMap(
            particle.mLatestMap, latestPoses, latestScanData);
    }

    /* Inspect the memory usage for the latest grid maps */
    this->mMetrics.mLatestMapMemoryUsage->Observe(
        this->InspectLatestMapMemoryUsage());
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
            return accWeight + particle.mWeight; });
    const double weightUnit =
        weightSum / static_cast<double>(numOfParticles);

    /* Perform low-variance sampling */
    std::vector<int> resampledIndices;
    resampledIndices.reserve(numOfParticles);

    double accWeight = uniformDist(this->mRandEngine) * weightUnit;
    double accParticleWeight = this->mParticles[0].mWeight;
    int particleIdx = 0;

    for (int i = 0; i < static_cast<int>(numOfParticles); ++i) {
        while (accWeight > accParticleWeight)
            accParticleWeight += this->mParticles[++particleIdx].mWeight;

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
            return accWeight + particle.mWeight * particle.mWeight; });
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

/* Inspect the memory usage for the particle grid maps */
std::uint64_t GridMapBuilder::InspectGridMapMemoryUsage() const
{
    return std::accumulate(
        this->mParticles.begin(), this->mParticles.end(), 0,
        [](const std::uint64_t value, const Particle& particle) {
            return value + particle.mMap.InspectMemoryUsage(); });
}

/* Inspect the memory usage for the latest grid maps */
std::uint64_t GridMapBuilder::InspectLatestMapMemoryUsage() const
{
    return std::accumulate(
        this->mParticles.begin(), this->mParticles.end(), 0,
        [](const std::uint64_t value, const Particle& particle) {
            return value + particle.mLatestMap.InspectMemoryUsage(); });
}

/* Inspect the memory usage for the trajectory nodes */
std::uint64_t GridMapBuilder::InspectTrajectoryMemoryUsage() const
{
    std::uint64_t memoryUsage = 0;
    std::unordered_set<std::uintptr_t> nodeSet;

    for (std::size_t i = 0; i < this->mParticles.size(); ++i) {
        auto nodePtr = this->mParticles[i].mNode;

        while (nodePtr != nullptr) {
            auto nodeId = reinterpret_cast<std::uintptr_t>(nodePtr.get());

            if (nodeSet.find(nodeId) != nodeSet.end())
                break;

            memoryUsage += sizeof(nodePtr->StampedPose()) +
                           sizeof(nodePtr->Parent());
            nodeSet.insert(nodeId);
            nodePtr = nodePtr->Parent();
        }
    }

    return memoryUsage;
}

/*
 * Utility function implementations
 */

/* Compute the maximum of a 'winSize' pixel wide row at each pixel */
void SlidingWindowMaxRow(const GridMap& gridMap,
                         const int winSize,
                         const Point2D<int>& idxMin,
                         ConstMap& intermediateMap)
{
    /* Each grid cell in the grid map stores an occupancy probability value
     * which is discretized to the unsigned integer (std::uint16_t) */
    using ValueType = GridMap::GridType::ValueType;

    /* Compute the maximum for each column */
    const ValueType unknownValue = gridMap.UnknownValue();
    int colIdx = 0;

    std::function<ValueType(int)> inFunc =
        [&colIdx, &gridMap, &idxMin, unknownValue](int rowIdx) {
        return gridMap.ValueOr(idxMin.mY + rowIdx,
                               idxMin.mX + colIdx,
                               unknownValue); };

    std::function<void(int, ValueType)> outFunc =
        [&colIdx, &idxMin, &intermediateMap, unknownValue](
            int rowIdx, ValueType maxValue) {
        intermediateMap.SetValueUnchecked(rowIdx, colIdx, maxValue); };

    const int rows = intermediateMap.Rows();
    const int cols = intermediateMap.Cols();

    /* Apply the sliding window maximum function */
    for (colIdx = 0; colIdx < cols; ++colIdx)
        SlidingWindowMax(inFunc, outFunc, rows, winSize);
}

/* Compute the maximum of a 'winSize' pixel wide column at each pixel */
void SlidingWindowMaxCol(const ConstMap& intermediateMap,
                         const int winSize,
                         ConstMap& precompMap)
{
    /* Make sure that the resulting map has the same size
     * as the intermediate map */
    Assert(precompMap.BlockRows() == intermediateMap.BlockRows());
    Assert(precompMap.BlockCols() == intermediateMap.BlockCols());

    /* Each grid cell in the grid map stores an occupancy probability value
     * which is discretized to the unsigned integer (std::uint16_t) */
    using ValueType = GridMap::GridType::ValueType;

    /* Compute the maximum for each row */
    const ValueType unknownValue = intermediateMap.UnknownValue();
    int rowIdx = 0;

    std::function<ValueType(int)> inFunc =
        [&rowIdx, &intermediateMap, unknownValue](int colIdx) {
        return intermediateMap.ValueOr(rowIdx, colIdx, unknownValue); };

    std::function<void(int, ValueType)> outFunc =
        [&rowIdx, &precompMap, unknownValue](
            int colIdx, ValueType maxValue) {
        precompMap.SetValueUnchecked(rowIdx, colIdx, maxValue); };

    const int rows = precompMap.Rows();
    const int cols = precompMap.Cols();

    /* Apply the sliding window maximum function */
    for (rowIdx = 0; rowIdx < rows; ++rowIdx)
        SlidingWindowMax(inFunc, outFunc, cols, winSize);
}

/* Precompute grid map for efficiency */
void PrecomputeGridMap(const GridMap& gridMap,
                       const int winSize,
                       const Point2D<int>& idxMin,
                       ConstMap& intermediateMap,
                       ConstMap& precompMap)
{
    /* Check if the resulting grid map size is same as the intermediate map */
    Assert(intermediateMap.Rows() == precompMap.Rows());
    Assert(intermediateMap.Cols() == precompMap.Cols());
    /* Check if the minimum grid cell index is inside the grid map */
    Assert(idxMin.mY >= 0 && idxMin.mY < gridMap.Rows());
    Assert(idxMin.mX >= 0 && idxMin.mX < gridMap.Cols());
    /* Check if the maximum grid cell index is inside the grid map */
    Assert(idxMin.mY + precompMap.Rows() > 0);
    Assert(idxMin.mX + precompMap.Cols() > 0);

    /* Each pixel in the `precompMap` stores the maximum of the
     * occupancy probability values of 'winSize' * 'winSize'
     * box of pixels beginning there */

    /* Store the maximum of the 'winSize' pixel wide row */
    SlidingWindowMaxRow(gridMap, winSize, idxMin, intermediateMap);
    /* Store the maximum of the 'winSize' pixel wide column */
    SlidingWindowMaxCol(intermediateMap, winSize, precompMap);

    /* Adjust the coordinate systems to make `minPos` correspond to
     * the origin grid cell (0, 0) */
    const auto minPos = gridMap.IndexToPosition(idxMin.mY, idxMin.mX);
    intermediateMap.SetPosOffset(minPos.mX, minPos.mY);
    precompMap.SetPosOffset(minPos.mX, minPos.mY);
}

} /* namespace Mapping */
} /* namespace MyGMapping */
