
/* scan_matcher_hill_climbing.cpp */

#include "my_gmapping/mapping/scan_matcher_hill_climbing.hpp"

namespace MyGMapping {
namespace Mapping {

/* Constructor */
ScanMatcherHillClimbingMetrics::ScanMatcherHillClimbingMetrics(
    const std::string& scanMatcherName) :
    mOptimizationTime(nullptr),
    mDiffTranslation(nullptr),
    mDiffRotation(nullptr),
    mNumOfIterations(nullptr),
    mNumOfRefinements(nullptr),
    mLikelihoodValue(nullptr),
    mNumOfScans(nullptr)
{
    /* Retrieve the metrics manager instance */
    auto* const pMetricManager = Metric::MetricManager::Instance();

    /* Register the value sequence metrics */
    this->mOptimizationTime = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".OptimizationTime");

    this->mDiffTranslation = pMetricManager->AddValueSequence<float>(
        scanMatcherName + ".DiffTranslation");
    this->mDiffRotation = pMetricManager->AddValueSequence<float>(
        scanMatcherName + ".DiffRotation");
    this->mNumOfIterations = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".NumOfIterations");
    this->mNumOfRefinements = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".NumOfRefinements");

    this->mLikelihoodValue = pMetricManager->AddValueSequence<float>(
        scanMatcherName + ".CostValue");
    this->mNumOfScans = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".NumOfScans");
}

/* Reserve the buffer to store the processing times */
void ScanMatcherHillClimbingMetrics::Resize(
    const std::size_t numOfParticles)
{
    this->mTimes.resize(numOfParticles);
    this->mParams.resize(numOfParticles);
}

/* Set the processing times for each particle */
void ScanMatcherHillClimbingMetrics::SetTimes(
    const std::size_t idx, const Times& times)
{
    this->mTimes[idx] = times;
}

/* Set the parameter settings for each particle */
void ScanMatcherHillClimbingMetrics::SetParameters(
    const std::size_t idx, const Parameters& parameters)
{
    this->mParams[idx] = parameters;
}

/* Collect the particle-wise metrics and update the overall metrics */
void ScanMatcherHillClimbingMetrics::Update(
    const std::size_t bestIdx)
{
    auto collectTimes = [&](
        std::function<int(int, const Times&)> selector) {
        return std::accumulate(this->mTimes.begin(),
                               this->mTimes.end(), 0, selector); };

    /* Compute the sum of the processing times */
    this->mOptimizationTime->Observe(collectTimes(
        [](int value, const Times& times) {
            return value + times.mOptimizationTime; }));

    /* Store the metric values of the best particle */
    const auto& bestParticle = this->mParams[bestIdx];
    this->mDiffTranslation->Observe(bestParticle.mDiffTranslation);
    this->mDiffRotation->Observe(bestParticle.mDiffRotation);
    this->mNumOfIterations->Observe(bestParticle.mNumOfIterations);
    this->mNumOfRefinements->Observe(bestParticle.mNumOfRefinements);
    this->mLikelihoodValue->Observe(bestParticle.mLikelihoodValue);
    this->mNumOfScans->Observe(bestParticle.mNumOfScans);
}

/* Constructor */
ScanMatcherHillClimbing::ScanMatcherHillClimbing(
    const std::string& scanMatcherName,
    LikelihoodFuncPtr&& likelihoodFunc,
    const double linearDelta,
    const double angularDelta,
    const int maxIterations,
    const int numOfRefinements) :
    ScanMatcher(scanMatcherName),
    mLikelihoodFunc(std::move(likelihoodFunc)),
    mLinearDelta(linearDelta),
    mAngularDelta(angularDelta),
    mMaxIterations(maxIterations),
    mNumOfRefinements(numOfRefinements),
    mMetrics(scanMatcherName)
{
}

/* Optimize pose by scan matching methods */
ScanMatchingResult ScanMatcherHillClimbing::OptimizePoseCore(
    const std::size_t particleIdx,
    const ScanMatchingQuery& query,
    const Sensor::ScanDataPtr<double>& scanData)
{
    /* Create the timer */
    Metric::Timer timer;
    /* Metric values for each particle */
    ScanMatcherMetrics::Times times;
    ScanMatcherMetrics::Parameters params;

    const auto& gridMap = query.mGridMap;
    const auto& initialPose = query.mInitialPose;

    /* Constants used for gradient descent method */
    static const double moveX[] = { 1.0, -1.0, 0.0, 0.0, 0.0, 0.0 };
    static const double moveY[] = { 0.0, 0.0, 1.0, -1.0, 0.0, 0.0 };
    static const double moveTheta[] = { 0.0, 0.0, 0.0, 0.0, 1.0, -1.0 };

    /* Calculate the sensor pose from the initial robot pose */
    const RobotPose2D<double> relPose = scanData->RelativeSensorPose();
    const RobotPose2D<double> sensorPose = Compound(initialPose, relPose);

    /* Best score and best pose */
    double bestLikelihood = this->mLikelihoodFunc->Likelihood(
        gridMap, scanData, sensorPose);
    RobotPose2D<double> bestPose = sensorPose;

    /* The number of the refinements */
    int numOfIterations = 0;
    int numOfRefinements = 0;
    double currentLinearDelta = this->mLinearDelta;
    double currentAngularDelta = this->mAngularDelta;
    bool poseUpdated = false;

    do {
        /* Local optimization */
        double bestLocalLikelihood = bestLikelihood;
        RobotPose2D<double> bestLocalPose = bestPose;
        poseUpdated = false;

        for (int i = 0; i < 6; ++i) {
            /* Move forward, backward, left, right,
             * rotate left and rotate right a little bit
             * then calculate the scan matching score value */
            const RobotPose2D<double> localPose {
                bestPose.mX + moveX[i] * currentLinearDelta,
                bestPose.mY + moveY[i] * currentLinearDelta,
                bestPose.mTheta + moveTheta[i] * currentAngularDelta };

            /* Calculate the scan matching score */
            const double localLikelihood = this->mLikelihoodFunc->Likelihood(
                gridMap, scanData, localPose);

            if (bestLocalLikelihood < localLikelihood) {
                bestLocalLikelihood = localLikelihood;
                bestLocalPose = localPose;
                poseUpdated = true;
            }
        }

        /* Update best pose */
        if (poseUpdated) {
            bestLikelihood = bestLocalLikelihood;
            bestPose = bestLocalPose;
        } else {
            /* Update the delta value if pose not improved */
            ++numOfRefinements;
            currentLinearDelta *= 0.5;
            currentAngularDelta *= 0.5;
        }
    } while ((poseUpdated || numOfRefinements < this->mNumOfRefinements) &&
             (++numOfIterations < this->mMaxIterations));

    /* Compute the robot pose from the sensor pose */
    const RobotPose2D<double> estimatedPose = MoveBackward(bestPose, relPose);
    const RobotPose2D<double> diffPose =
        InverseCompound(initialPose, estimatedPose);

    /* Set the resulting likelihood */
    const double likelihood = bestLikelihood;
    const double normalizedLikelihood = likelihood / scanData->NumOfScans();
    /* Set the resulting score (same as the likelihood) */
    const double score = likelihood;
    const double normalizedScore = normalizedLikelihood;

    /* Update the metrics */
    times.mOptimizationTime = timer.ElapsedMicro();
    timer.Stop();

    params.mDiffTranslation = Distance(diffPose);
    params.mDiffRotation = std::abs(diffPose.mTheta);
    params.mNumOfIterations = numOfIterations;
    params.mNumOfRefinements = numOfRefinements;
    params.mLikelihoodValue = normalizedLikelihood;
    params.mNumOfScans = scanData->NumOfScans();

    /* Set the metrics */
    this->mMetrics.SetTimes(particleIdx, times);
    this->mMetrics.SetParameters(particleIdx, params);

    return ScanMatchingResult { initialPose, estimatedPose,
                                normalizedLikelihood, likelihood,
                                normalizedScore, score };
}

/* Optimize poses for multiple particles by scan matching methods */
ScanMatchingResultVector ScanMatcherHillClimbing::OptimizePose(
    const ScanMatchingQueryVector& queries,
    const Sensor::ScanDataPtr<double>& scanData)
{
    ScanMatchingResultVector results;
    results.resize(queries.size());

    /* Resize the buffer to store the metrics for each particle */
    this->mMetrics.Resize(queries.size());

    /* Optimize the pose for each particle (parallelized using OpenMP) */
#pragma omp parallel for
    for (std::size_t i = 0; i < queries.size(); ++i)
        results[i] = this->OptimizePoseCore(i, queries[i], scanData);

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

} /* namespace Mapping */
} /* namespace MyGMapping */
