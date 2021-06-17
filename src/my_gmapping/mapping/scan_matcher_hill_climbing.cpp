
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

    /* Register the distribution metrics */
    this->mOptimizationTime = pMetricManager->AddDistribution(
        scanMatcherName + ".OptimizationTime");
    this->mDiffTranslation = pMetricManager->AddDistribution(
        scanMatcherName + ".DiffTranslation");
    this->mDiffRotation = pMetricManager->AddDistribution(
        scanMatcherName + ".DiffRotation");
    this->mNumOfIterations = pMetricManager->AddDistribution(
        scanMatcherName + ".NumOfIterations");
    this->mNumOfRefinements = pMetricManager->AddDistribution(
        scanMatcherName + ".NumOfRefinements");

    /* Register the value sequence metrics */
    this->mLikelihoodValue = pMetricManager->AddValueSequence<float>(
        scanMatcherName + ".CostValue");
    this->mNumOfScans = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".NumOfScans");
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
    const ScanMatchingQuery& query,
    const Sensor::ScanDataPtr<double>& scanData)
{
    /* Create the timer */
    Metric::Timer timer;

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

    /* Set the resulting likelihood */
    const double likelihood = bestLikelihood;
    const double normalizedLikelihood = likelihood / scanData->NumOfScans();
    /* Set the resulting score (same as the likelihood) */
    const double score = likelihood;
    const double normalizedScore = normalizedLikelihood;

    /* Update the metrics */
    this->mMetrics.mOptimizationTime->Observe(timer.ElapsedMicro());
    this->mMetrics.mDiffTranslation->Observe(
        Distance(initialPose, estimatedPose));
    this->mMetrics.mDiffRotation->Observe(
        std::abs(initialPose.mTheta - estimatedPose.mTheta));
    this->mMetrics.mNumOfIterations->Observe(numOfIterations);
    this->mMetrics.mNumOfRefinements->Observe(numOfRefinements);

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

    /* Optimize the pose for each particle (parallelized using OpenMP) */
#pragma omp parallel for
    for (std::size_t i = 0; i < queries.size(); ++i)
        results[i] = this->OptimizePoseCore(queries[i], scanData);

    /* Determine the maximum likelihood value and its corresponding score */
    const auto bestIt = std::max_element(
        results.begin(), results.end(),
        [](const ScanMatchingResult& lhs, const ScanMatchingResult& rhs) {
            return lhs.mLikelihood < rhs.mLikelihood; });

    /* Update the normalized likelihood value of the best particle */
    this->mMetrics.mLikelihoodValue->Observe(bestIt->mNormalizedLikelihood);
    /* Update the number of the scan points */
    this->mMetrics.mNumOfScans->Observe(scanData->NumOfScans());

    return results;
}

} /* namespace Mapping */
} /* namespace MyGMapping */
