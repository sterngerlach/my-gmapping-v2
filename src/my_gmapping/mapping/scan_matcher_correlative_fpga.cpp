
/* scan_matcher_correlative_fpga.cpp */

#include "my_gmapping/mapping/scan_matcher_correlative_fpga.hpp"

namespace MyGMapping {
namespace Mapping {

using namespace Hardware;

/* Static assertions for safety */
static_assert(std::is_same<GridMapType::StorageType, std::uint16_t>::value,
              "Grid map should store occupancy probability values as "
              "discretized 16-bit unsigned integers");
static_assert(GridMapType::CellType::UnknownRaw == 0,
              "Unknown occupancy probability should be defined as 0");

/* Constructor */
ScanMatcherCorrelativeFPGAMetrics::ScanMatcherCorrelativeFPGAMetrics(
    const std::string& scanMatcherName) :
    mInputSetupTime(nullptr),
    mSetupIPTime(nullptr),
    mScanSendTime(nullptr),
    mMapSendTime(nullptr),
    mOptimizationTime(nullptr),
    mWaitIPTime(nullptr),
    mScanMatchingTime(nullptr),
    mDiffTranslation(nullptr),
    mDiffRotation(nullptr),
    mWinSizeX(nullptr),
    mWinSizeY(nullptr),
    mWinSizeTheta(nullptr),
    mStepSizeX(nullptr),
    mStepSizeY(nullptr),
    mStepSizeTheta(nullptr),
    mMapSizeX(nullptr),
    mMapSizeY(nullptr),
    mMapChunks(nullptr),
    mScanTransferSkip(nullptr),
    mMapTransferSkip(nullptr),
    mScoreValue(nullptr),
    mLikelihoodValue(nullptr),
    mNumOfTransferredScans(nullptr)
{
    /* Retrieve the metrics manager instance */
    auto* const pMetricManager = Metric::MetricManager::Instance();

    /* Register the counter metrics */
    this->mMapChunks = pMetricManager->AddCounter(
        scanMatcherName + ".MapChunks");
    this->mScanTransferSkip = pMetricManager->AddCounter(
        scanMatcherName + ".ScanTransferSkip");
    this->mMapTransferSkip = pMetricManager->AddCounter(
        scanMatcherName + ".MapTransferSkip");

    /* Register the distribution metrics */
    this->mInputSetupTime = pMetricManager->AddDistribution(
        scanMatcherName + ".InputSetupTime");
    this->mSetupIPTime = pMetricManager->AddDistribution(
        scanMatcherName + ".SetupIPTime");
    this->mScanSendTime = pMetricManager->AddDistribution(
        scanMatcherName + ".ScanSendTime");
    this->mMapSendTime = pMetricManager->AddDistribution(
        scanMatcherName + ".MapSendTime");
    this->mOptimizationTime = pMetricManager->AddDistribution(
        scanMatcherName + ".OptimizationTime");
    this->mWaitIPTime = pMetricManager->AddDistribution(
        scanMatcherName + ".WaitIPTime");
    this->mScanMatchingTime = pMetricManager->AddDistribution(
        scanMatcherName + ".ScanMatchingTime");
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

    /* Register the histogram metrics */
    const Metric::BucketBoundaries mapSizeBuckets =
        Metric::Histogram::CreateFixedWidthBoundaries(0.0, 400.0, 20.0);
    this->mMapSizeX = pMetricManager->AddHistogram(
        scanMatcherName + ".MapSizeX", mapSizeBuckets);
    this->mMapSizeY = pMetricManager->AddHistogram(
        scanMatcherName + ".MapSizeY", mapSizeBuckets);

    /* Register the value sequence metrics */
    this->mScoreValue = pMetricManager->AddValueSequence<float>(
        scanMatcherName + ".ScoreValue");
    this->mLikelihoodValue = pMetricManager->AddValueSequence<float>(
        scanMatcherName + ".LikelihoodValue");
    this->mNumOfTransferredScans = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".NumOfTransferredScans");
}

/* Constructor */
ScanMatcherCorrelativeFPGA::ScanMatcherCorrelativeFPGA(
    const std::string& scanMatcherName,
    LikelihoodFuncPtr&& likelihoodFunc,
    const double rangeX,
    const double rangeY,
    const double rangeTheta,
    const double scanRangeMax,
    IPConfig&& scanMatcherConfig0,
    IPConfig&& scanMatcherConfig1,
    IPCommonConfig&& commonConfig,
    AxiDmaConfig&& axiDmaConfig0,
    AxiDmaConfig&& axiDmaConfig1) :
    ScanMatcher(scanMatcherName),
    mLikelihoodFunc(std::move(likelihoodFunc)),
    mRangeX(rangeX),
    mRangeY(rangeY),
    mRangeTheta(rangeTheta),
    mScanRangeMax(scanRangeMax),
    mConfig({ std::move(scanMatcherConfig0),
              std::move(scanMatcherConfig1) }),
    mCommonConfig(std::move(commonConfig)),
    mAxiDmaConfig({ std::move(axiDmaConfig0),
                    std::move(axiDmaConfig1) }),
    mMetrics(scanMatcherName)
{
    /* Check the information of the scan matcher IP core */
    XAssert(this->mCommonConfig.mMaxNumOfScans > 0,
            "Maximum number of the scan points should be greater than 0");
    XAssert(this->mCommonConfig.mMapResolution > 0.0,
            "Grid map resolution should be greater than 0");
    XAssert(this->mCommonConfig.mMaxMapSizeX > 0,
            "Maximum width of the grid map (in the number of grid cells) "
            "should be greater than 0");
    XAssert(this->mCommonConfig.mMaxMapSizeY > 0,
            "Maximum height of the grid map (in the number of grid cells) "
            "should be greater than 0");
    XAssert(this->mCommonConfig.mLowResolution > 1,
            "Resolution of the coarse grid map (in the number of grid cells) "
            "should be greater than 0");
    XAssert(this->mCommonConfig.mMapBitWidth > 0 &&
            this->mCommonConfig.mMapBitWidth <= 8,
            "Bit width of the discretized occupancy probability "
            "should be in the range between 1 and 8");
    XAssert(this->mCommonConfig.mMapChunkWidth > 0 &&
            this->mCommonConfig.mMapChunkWidth <= 8,
            "Width of the grid map chunk (consecutive grid map cells) "
            "should be in the range between 1 and 8");

    /* Initialize the scan matcher IP cores */
    this->Initialize(0);
    this->Initialize(1);
}

/* Optimize the particle poses based on the correlative scan matching */
void ScanMatcherCorrelativeFPGA::OptimizePose(
    const std::size_t numOfParticles,
    const std::vector<const GridMapType*>& particleMaps,
    const Sensor::ScanDataPtr<double>& scanData,
    const std::vector<RobotPose2D<double>>& initialPoses,
    std::vector<RobotPose2D<double>>& estimatedPoses,
    std::vector<double>& likelihoodValues)
{
    /* Check the number of the particles */
    XAssert(numOfParticles % NumOfIPCores == 0,
            "Number of the particles must be multiple of the "
            "number of the scan matcher IP cores implemented on the device");
    XAssert(particleMaps.size() == numOfParticles,
            "Mismatch between the number of the grid maps and "
            "the number of the particles");
    XAssert(initialPoses.size() == numOfParticles,
            "Mismatch between the number of the initial poses and "
            "the number of the particles");

    /* Reserve the vector to store the results */
    estimatedPoses.clear();
    estimatedPoses.reserve(numOfParticles);
    likelihoodValues.clear();
    likelihoodValues.reserve(numOfParticles);

    std::vector<double> normalizedLikelihoods;
    normalizedLikelihoods.reserve(numOfParticles);
    std::vector<double> normalizedScores;
    normalizedScores.reserve(numOfParticles);

    const std::size_t halfNumOfParticles = numOfParticles / 2;

    std::vector<RobotPose2D<double>> estimatedPoses1;
    estimatedPoses1.reserve(halfNumOfParticles);
    std::vector<double> likelihoodValues1;
    likelihoodValues1.reserve(halfNumOfParticles);
    std::vector<double> normalizedLikelihoods1;
    normalizedLikelihoods1.reserve(halfNumOfParticles);
    std::vector<double> normalizedScores1;
    normalizedScores1.reserve(halfNumOfParticles);

    /* Start a sub-thread to process the second half of the particles */
    this->mSubThread = std::thread([&]() {
        this->OptimizePoseCore(1, halfNumOfParticles, numOfParticles,
                               particleMaps, scanData, initialPoses,
                               estimatedPoses1, likelihoodValues1,
                               normalizedLikelihoods1, normalizedScores1); });

    /* Start to process the first half of the particles */
    this->OptimizePoseCore(0, 0, halfNumOfParticles,
                           particleMaps, scanData, initialPoses,
                           estimatedPoses, likelihoodValues,
                           normalizedLikelihoods, normalizedScores);

    /* Wait for a sub-thread to finish */
    if (this->mSubThread.joinable())
        this->mSubThread.join();

    /* Concatenate the results */
    estimatedPoses.insert(
        estimatedPoses.end(),
        std::make_move_iterator(estimatedPoses1.begin()),
        std::make_move_iterator(estimatedPoses1.end()));
    likelihoodValues.insert(
        likelihoodValues.end(),
        std::make_move_iterator(likelihoodValues1.begin()),
        std::make_move_iterator(likelihoodValues1.end()));
    normalizedLikelihoods.insert(
        normalizedLikelihoods.end(),
        std::make_move_iterator(normalizedLikelihoods1.begin()),
        std::make_move_iterator(normalizedLikelihoods1.end()));
    normalizedScores.insert(
        normalizedScores.end(),
        std::make_move_iterator(normalizedScores1.begin()),
        std::make_move_iterator(normalizedScores1.end()));

    /* Determine the maximum likelihood value and its corresponding score */
    const auto bestIt = std::max_element(
        likelihoodValues.begin(), likelihoodValues.end());
    const auto bestIdx = std::distance(likelihoodValues.begin(), bestIt);

    /* Determine the number of the scan points to be transferred */
    const int maxNumOfScans = this->mCommonConfig.mMaxNumOfScans;
    const int numOfScans = static_cast<int>(scanData->NumOfScans());
    const int numOfScansTransferred = std::min(maxNumOfScans, numOfScans);

    /* Update the metrics */
    this->mMetrics.mScoreValue->Observe(normalizedScores[bestIdx]);
    this->mMetrics.mLikelihoodValue->Observe(normalizedLikelihoods[bestIdx]);
    this->mMetrics.mNumOfTransferredScans->Observe(numOfScansTransferred);
}

/* Optimize the particle poses using the scan matcher IP core */
void ScanMatcherCorrelativeFPGA::OptimizePoseCore(
    const int coreId,
    const std::size_t particleIdxBegin,
    const std::size_t particleIdxEnd,
    const std::vector<const GridMapType*>& particleMaps,
    const Sensor::ScanDataPtr<double>& scanData,
    const std::vector<RobotPose2D<double>>& initialPoses,
    std::vector<RobotPose2D<double>>& estimatedPoses,
    std::vector<double>& likelihoodValues,
    std::vector<double>& normalizedLikelihoodValues,
    std::vector<double>& normalizedScores)
{
    /* Determine the number of the scan points to be transferred */
    const int maxNumOfScans = this->mCommonConfig.mMaxNumOfScans;
    const int numOfScans = static_cast<int>(scanData->NumOfScans());
    const int numOfScansTransferred = std::min(maxNumOfScans, numOfScans);

    /* Determine the maximum scan range considered valid */
    const float scanRangeMax = static_cast<float>(
        std::min(this->mScanRangeMax, scanData->MaxRange()));
    /* Score threshold is always set to zero */
    const int quantizedScoreThreshold = 0;

    /* Determine the search step */
    double stepX;
    double stepY;
    double stepTheta;
    this->ComputeSearchStep(this->mCommonConfig.mMapResolution, scanData,
                            stepX, stepY, stepTheta);

    /* Determine the radius of the search window */
    const int originalWinX = static_cast<int>(
        std::ceil(0.5 * this->mRangeX / stepX));
    const int originalWinY = static_cast<int>(
        std::ceil(0.5 * this->mRangeY / stepY));
    const int originalWinTheta = static_cast<int>(
        std::ceil(0.5 * this->mRangeTheta / stepTheta));

    /* Determine the actual radius of the search window */
    const int winX = std::clamp(originalWinX, 1,
                                this->mCommonConfig.mMaxMapSizeX / 2);
    const int winY = std::clamp(originalWinY, 1,
                                this->mCommonConfig.mMaxMapSizeY / 2);
    const int winTheta = std::max(1, originalWinTheta);

    /* Optimize the particle poses for each particle */
    for (std::size_t i = particleIdxBegin; i < particleIdxEnd; ++i) {
        /* Create the timer */
        Metric::Timer outerTimer;
        Metric::Timer timer;

        /* Compute the sensor pose from the initial particle pose */
        const RobotPose2D<double> sensorPose =
            Compound(initialPoses[i], scanData->RelativeSensorPose());

        /* Compute the minimum possible sensor pose */
        const RobotPose2D<double> minSensorPose {
            sensorPose.mX - stepX * winX,
            sensorPose.mY - stepY * winY,
            sensorPose.mTheta - stepTheta * winTheta };

        /* Desired size of the grid map */
        const int maxMapSizeX = this->mCommonConfig.mMaxMapSizeX;
        const int maxMapSizeY = this->mCommonConfig.mMaxMapSizeY;

        /* Compute the center position of the cropped grid map */
        const Point2D<int> gridMapCenterIdx =
            particleMaps[i]->MapCoordinateToCellIndex(
                sensorPose.mX, sensorPose.mY);
        /* Use std::clamp() here since the above `gridMapCenterIdx` could be
         * out-of-bounds (the current particle is outside of the grid map) */
        const Point2D<int> minimumPossibleIdx {
            std::clamp(gridMapCenterIdx.mX - maxMapSizeX,
                       0, particleMaps[i]->NumCellsX() - 1),
            std::clamp(gridMapCenterIdx.mY - maxMapSizeY,
                       0, particleMaps[i]->NumCellsY() - 1) };
        const Point2D<int> maximumPossibleIdx {
            std::clamp(gridMapCenterIdx.mX + maxMapSizeX,
                       1, particleMaps[i]->NumCellsX()),
            std::clamp(gridMapCenterIdx.mY + maxMapSizeY,
                       1, particleMaps[i]->NumCellsY()) };
        const Point2D<int> finalCenterIdx {
            (maximumPossibleIdx.mX + minimumPossibleIdx.mX) / 2,
            (maximumPossibleIdx.mY + minimumPossibleIdx.mY) / 2 };

        /* Crop the grid map around the center position */
        const Point2D<int> gridMapMinIdx {
            std::max(finalCenterIdx.mX - maxMapSizeX / 2, 0),
            std::max(finalCenterIdx.mY - maxMapSizeY / 2, 0) };
        const Point2D<int> gridMapMaxIdx {
            std::min(finalCenterIdx.mX + maxMapSizeX / 2,
                     particleMaps[i]->NumCellsX()),
            std::min(finalCenterIdx.mY + maxMapSizeY / 2,
                     particleMaps[i]->NumCellsY()) };

        /* Compute the size of the cropped grid map */
        const Point2D<int> gridMapSize {
            gridMapMaxIdx.mX - gridMapMinIdx.mX,
            gridMapMaxIdx.mY - gridMapMinIdx.mY };
        /* Compute the minimum position of the cropped grid map */
        const Point2D<double> gridMapMinPos =
            particleMaps[i]->CellIndexToMapCoordinate(gridMapMinIdx);

        /* Update the metrics and restart the timer */
        this->mMetrics.mInputSetupTime->Observe(timer.ElapsedMicro());
        timer.Start();

        /* Set the scan matching parameters through AXI4-Lite slave interface */
        this->SetParameterRegisters(
            coreId, numOfScansTransferred,
            scanRangeMax, quantizedScoreThreshold,
            minSensorPose, gridMapSize, gridMapMinPos,
            winX * 2, winY * 2, winTheta * 2, stepX, stepY, stepTheta);
        /* Start the scan matcher IP core */
        this->StartIPCore(coreId);

        /* Update the metrics and restart the timer */
        this->mMetrics.mSetupIPTime->Observe(timer.ElapsedMicro());
        timer.Start();

        /* Send the scan data for the first particle only */
        const bool scanDataTransferred = (i == particleIdxBegin);
        this->SendScanData(coreId, scanDataTransferred, scanData);
        /* Update the metrics and restart the timer */
        this->mMetrics.mScanSendTime->Observe(timer.ElapsedMicro());
        timer.Start();

        /* Send the grid map */
        this->SendGridMap(coreId, *particleMaps[i], gridMapMinIdx, gridMapSize);
        /* Update the metrics and restart the timer */
        this->mMetrics.mMapSendTime->Observe(timer.ElapsedMicro());
        timer.Start();

        /* Receive the result */
        int scoreMax;
        int bestX;
        int bestY;
        int bestTheta;
        this->ReceiveResult(coreId, scoreMax, bestX, bestY, bestTheta);
        /* Update the metrics and restart the timer */
        this->mMetrics.mOptimizationTime->Observe(timer.ElapsedMicro());
        timer.Start();

        /* Wait for the scan matcher IP core */
        this->WaitIPCore(coreId);
        /* Update the metrics and stop the timer */
        this->mMetrics.mWaitIPTime->Observe(timer.ElapsedMicro());

        /* The appropriate solution is found if the maximum score is
         * larger than (not larger than or equal to) the score threshold */
        const bool poseFound = scoreMax > quantizedScoreThreshold;
        const double normalizedScore =
            static_cast<double>(scoreMax) /
            static_cast<double>(numOfScansTransferred) /
            static_cast<double>((1 << this->mCommonConfig.mMapBitWidth) - 1);
        /* Compute the best sensor pose */
        const RobotPose2D<double> bestSensorPose {
            minSensorPose.mX + bestX * stepX,
            minSensorPose.mY + bestY * stepY,
            minSensorPose.mTheta + bestTheta * stepTheta };

        /* Compute the estimated robot pose and the likelihood value */
        RobotPose2D<double> estimatedPose =
            MoveBackward(bestSensorPose, scanData->RelativeSensorPose());
        const double likelihoodValue = this->mLikelihoodFunc->Likelihood(
            *particleMaps[i], scanData, bestSensorPose);
        const double normalizedLikelihood =
            likelihoodValue / scanData->NumOfScans();

        /* Update the metrics */
        this->mMetrics.mScanMatchingTime->Observe(outerTimer.ElapsedMicro());
        this->mMetrics.mDiffTranslation->Observe(
            Distance(initialPoses[i], estimatedPose));
        this->mMetrics.mDiffRotation->Observe(
            std::abs(initialPoses[i].mTheta - estimatedPose.mTheta));
        this->mMetrics.mMapSizeX->Observe(gridMapSize.mX);
        this->mMetrics.mMapSizeY->Observe(gridMapSize.mY);

        estimatedPoses.emplace_back(std::move(estimatedPose));
        likelihoodValues.emplace_back(likelihoodValue);
        normalizedLikelihoodValues.emplace_back(normalizedLikelihood);
        normalizedScores.emplace_back(normalizedScore);
    }

    /* Update the metrics */
    this->mMetrics.mWinSizeX->Observe(winX);
    this->mMetrics.mWinSizeY->Observe(winY);
    this->mMetrics.mWinSizeTheta->Observe(winTheta);
    this->mMetrics.mStepSizeX->Observe(stepX);
    this->mMetrics.mStepSizeY->Observe(stepY);
    this->mMetrics.mStepSizeTheta->Observe(stepTheta);
}

/* Compute the search step */
void ScanMatcherCorrelativeFPGA::ComputeSearchStep(
    const double mapResolution,
    const Sensor::ScanDataPtr<double>& scanData,
    double& stepX,
    double& stepY,
    double& stepTheta) const
{
    /* Determine the search step */
    const auto maxRangeIt = std::max_element(
        scanData->Ranges().cbegin(), scanData->Ranges().cend());
    const double maxRange = std::min(*maxRangeIt, this->mScanRangeMax);
    const double theta = mapResolution / maxRange;

    stepX = mapResolution;
    stepY = mapResolution;
    stepTheta = std::acos(1.0 - 0.5 * theta * theta);

    return;
}

/* Initialize the scan matcher IP core */
void ScanMatcherCorrelativeFPGA::Initialize(const int coreId)
{
    /* Ensure that the base address and address range are aligned */
    XAssert(this->mConfig[coreId].mAxiLiteSBaseAddress % 4 == 0,
            "Base address of the scan matcher IP core control registers "
            "must be 4-byte aligned");
    XAssert(this->mConfig[coreId].mAxiLiteSAddressRange % 4 == 0,
            "Address range of the scan matcher IP core control registers "
            "must be multiple of 4");

    /* Initialize the scan matcher IP core */
    this->mControlRegisters[coreId] = std::make_unique<MemoryMappedIO>(
        this->mConfig[coreId].mAxiLiteSBaseAddress,
        this->mConfig[coreId].mAxiLiteSAddressRange);

    /* Initialize the AXI DMA IP core */
    SharedMemoryMappedIO axiDmaRegisters {
        this->mAxiDmaConfig[coreId].mBaseAddress,
        this->mAxiDmaConfig[coreId].mAddressRange };
    this->mAxiDma[coreId] = std::make_unique<AxiSimpleDMA>(axiDmaRegisters);

    /* Initialize the CMA memory to store the input data */
    this->InitializeCMAMemoryInput(coreId);
    /* Initialize the CMA memory to store the output data */
    this->InitializeCMAMemoryOutput(coreId);

    /* Reset and halt the AXI DMA IP core */
    this->mAxiDma[coreId]->SendChannel().Reset();
    this->mAxiDma[coreId]->RecvChannel().Reset();
    this->mAxiDma[coreId]->SendChannel().Stop();
    this->mAxiDma[coreId]->RecvChannel().Stop();

    /* Start the AXI DMA transfer */
    this->mAxiDma[coreId]->SendChannel().Start();
    this->mAxiDma[coreId]->RecvChannel().Start();
}

/* Initialize the CMA memory for the input data */
void ScanMatcherCorrelativeFPGA::InitializeCMAMemoryInput(const int coreId)
{
    /* Compute the number of bytes required for transferring scan data
     * Each 64-bit element packs the scan range (32-bit floating point) and
     * the scan angle (32-bit floating point) */
    /* Add 1 to transfer the flag that indicates whether the scan data
     * should be transferred (scan data is only transferred for the first
     * particle and is reused for other particles) */
    const std::size_t scanDataLength = this->mCommonConfig.mMaxNumOfScans + 1;

    const int maxMapSizeX = this->mCommonConfig.mMaxMapSizeX;
    const int maxMapSizeY = this->mCommonConfig.mMaxMapSizeY;
    const int mapChunkWidth = this->mCommonConfig.mMapChunkWidth;

    /* Compute the number of bytes required for transferring grid map
     * Each 64-bit element (grid map chunk) packs 8 occupancy probability
     * values for 8 consecutive grid cells, which are quantized to
     * 8-bit unsigned integers when transferred */
    const std::size_t numOfMapChunksPerRow =
        (maxMapSizeX + mapChunkWidth - 1) / mapChunkWidth;
    /* Add 1 to transfer the flag that indicates whether the grid map
     * should be transferred (grid map is cached on the BRAM to reduce
     * the expensive data transfer cost) */
    const std::size_t gridMapLength = maxMapSizeY * numOfMapChunksPerRow + 1;

    /* The number of bytes required to store the input data */
    const std::size_t lengthInBytes =
        std::max(scanDataLength, gridMapLength) * sizeof(std::uint64_t);
    /* Allocate the CMA memory for the input data */
    this->mInputData[coreId].Initialize(
        static_cast<std::uint32_t>(lengthInBytes));
}

/* Initialize the CMA memory for the output data */
void ScanMatcherCorrelativeFPGA::InitializeCMAMemoryOutput(const int coreId)
{
    /* The number of bytes required to store the output data
     * The first 64-bit element contains the maximum score and the
     * discretized X-coordinate of the final sensor pose, and the
     * second 64-bit element contains the discretized Y-coordinate and
     * the discretized rotation angle of the final sensor pose */
    const std::size_t lengthInBytes = 2 * sizeof(std::uint64_t);
    /* Allocate the CMA memory for the output data */
    this->mOutputData[coreId].Initialize(
        static_cast<std::uint32_t>(lengthInBytes));
}

/* Set the scan matching parameters through AXI4-Lite slave interface */
void ScanMatcherCorrelativeFPGA::SetParameterRegisters(
    const int coreId,
    const int numOfScans,
    const double scanRangeMax,
    const int scoreThreshold,
    const RobotPose2D<double>& minSensorPose,
    const Point2D<int>& gridMapSize,
    const Point2D<double>& gridMapMinPos,
    const int winX, const int winY, const int winTheta,
    const double stepX, const double stepY, const double stepTheta)
{
    /* Set the actual number of the scan points */
    this->mControlRegisters[coreId]->Write(
        this->mCommonConfig.mAxiLiteSNumOfScans,
        Int32ToUInt32(numOfScans));
    /* Set the maximum scan range considered valid */
    this->mControlRegisters[coreId]->Write(
        this->mCommonConfig.mAxiLiteSScanRangeMax,
        FloatToUInt32(static_cast<float>(scanRangeMax)));
    /* Set the score threshold (for loop detection) */
    this->mControlRegisters[coreId]->Write(
        this->mCommonConfig.mAxiLiteSScoreThreshold,
        Int32ToUInt32(scoreThreshold));

    /* Set the minimum possible sensor pose */
    this->mControlRegisters[coreId]->Write(
        this->mCommonConfig.mAxiLiteSPoseX,
        FloatToUInt32(static_cast<float>(minSensorPose.mX)));
    this->mControlRegisters[coreId]->Write(
        this->mCommonConfig.mAxiLiteSPoseY,
        FloatToUInt32(static_cast<float>(minSensorPose.mY)));
    this->mControlRegisters[coreId]->Write(
        this->mCommonConfig.mAxiLiteSPoseTheta,
        FloatToUInt32(static_cast<float>(minSensorPose.mTheta)));

    /* Set the actual size of the cropped grid map */
    this->mControlRegisters[coreId]->Write(
        this->mCommonConfig.mAxiLiteSMapSizeX, Int32ToUInt32(gridMapSize.mX));
    this->mControlRegisters[coreId]->Write(
        this->mCommonConfig.mAxiLiteSMapSizeY, Int32ToUInt32(gridMapSize.mY));
    /* Set the minimum coordinate of the cropped grid map */
    this->mControlRegisters[coreId]->Write(
        this->mCommonConfig.mAxiLiteSMapMinX,
        FloatToUInt32(static_cast<float>(gridMapMinPos.mX)));
    this->mControlRegisters[coreId]->Write(
        this->mCommonConfig.mAxiLiteSMapMinY,
        FloatToUInt32(static_cast<float>(gridMapMinPos.mY)));

    /* Set the size of the search window */
    this->mControlRegisters[coreId]->Write(
        this->mCommonConfig.mAxiLiteSWinX, Int32ToUInt32(winX));
    this->mControlRegisters[coreId]->Write(
        this->mCommonConfig.mAxiLiteSWinY, Int32ToUInt32(winY));
    this->mControlRegisters[coreId]->Write(
        this->mCommonConfig.mAxiLiteSWinTheta, Int32ToUInt32(winTheta));
    /* Set the search step */
    this->mControlRegisters[coreId]->Write(
        this->mCommonConfig.mAxiLiteSStepX,
        FloatToUInt32(static_cast<float>(stepX)));
    this->mControlRegisters[coreId]->Write(
        this->mCommonConfig.mAxiLiteSStepY,
        FloatToUInt32(static_cast<float>(stepY)));
    this->mControlRegisters[coreId]->Write(
        this->mCommonConfig.mAxiLiteSStepTheta,
        FloatToUInt32(static_cast<float>(stepTheta)));
}

/* Send the scan data through AXI DMA */
void ScanMatcherCorrelativeFPGA::SendScanData(
    const int coreId,
    const bool scanDataTransferred,
    const Sensor::ScanDataPtr<double>& scanData)
{
    /* Retrieve the pointer to the CMA memory */
    volatile std::uint64_t* pInput =
        this->mInputData[coreId].Ptr<volatile std::uint64_t>();

    /* Do not transfer the scan data if it is redundant */
    if (!scanDataTransferred) {
        /* Write the flag to indicate that the scan data is not transferred */
        *pInput++ = static_cast<std::uint64_t>(0);

        /* Transfer the flag only using the AXI DMA IP core */
        this->mAxiDma[coreId]->SendChannel().Transfer(
            sizeof(std::uint64_t), this->mInputData[coreId].PhysicalAddress());

        /* Wait for the data transfer to complete if necessary */
        if (this->mCommonConfig.mWaitForDmaTransfer)
            this->mAxiDma[coreId]->SendChannel().Wait();

        /* Update the metrics */
        this->mMetrics.mScanTransferSkip->Increment();

        return;
    }

    /* Write the flag to indicate that the scan data is transferred */
    *pInput++ = static_cast<std::uint64_t>(1);

    const std::size_t maxNumOfScans =
        static_cast<std::size_t>(this->mCommonConfig.mMaxNumOfScans);
    const std::size_t numOfScans = scanData->NumOfScans();

    /* Subsample the scan points if the number of the scan points
     * exceeds the maximum that the scan matcher IP core can handle */
    if (maxNumOfScans < numOfScans) {
        /* Interval is greater than 1 and thus the scan points are sampled */
        const double idxInterval =
            static_cast<double>(numOfScans - 1) /
            static_cast<double>(maxNumOfScans - 1);
        double subsampledIdx = 0.0;

        for (std::size_t i = 0; i < maxNumOfScans; ++i) {
            /* Write the subsampled scan point */
            const std::size_t scanIdx =
                static_cast<std::size_t>(subsampledIdx);
            const float scanRange = static_cast<float>(
                scanData->RangeAt(scanIdx));
            const float scanAngle = static_cast<float>(
                scanData->AngleAt(scanIdx));
            *pInput++ = PackFloat(scanRange, scanAngle);
            /* Advance the index of the subsampled scan point */
            subsampledIdx += idxInterval;
        }
    } else {
        /* Write all the scan points */
        for (std::size_t i = 0; i < numOfScans; ++i) {
            /* Write the scan range and angle */
            const float scanRange = static_cast<float>(scanData->RangeAt(i));
            const float scanAngle = static_cast<float>(scanData->AngleAt(i));
            *pInput++ = PackFloat(scanRange, scanAngle);
        }
    }

    /* Transfer the scan points using the AXI DMA IP core */
    const std::size_t numOfScansTransferred =
        std::min(maxNumOfScans, numOfScans);
    const std::size_t transferLengthInBytes =
        (numOfScansTransferred + 1) * sizeof(std::uint64_t);
    this->mAxiDma[coreId]->SendChannel().Transfer(
        transferLengthInBytes, this->mInputData[coreId].PhysicalAddress());

    /* Wait for the data transfer to complete by polling the status register
     * of the AXI DMA IP core */
    if (this->mCommonConfig.mWaitForDmaTransfer)
        this->mAxiDma[coreId]->SendChannel().Wait();
}

/* Send the grid map through AXI DMA */
void ScanMatcherCorrelativeFPGA::SendGridMap(
    const int coreId,
    const GridMapType& gridMap,
    const Point2D<int>& gridMapMinIdx,
    const Point2D<int>& gridMapSize)
{
    /* Retrieve the pointer to the CMA memory */
    volatile std::uint64_t* pInput =
        this->mInputData[coreId].Ptr<volatile std::uint64_t>();

    /* Write the flag to indicate that the grid map is transferred */
    *pInput++ = static_cast<std::uint64_t>(1);

    /* Write the cropped grid map */
    const std::uint16_t unknownRawValue = gridMap.UnknownRawValue();
    const int mapChunkWidth = this->mCommonConfig.mMapChunkWidth;
    const int numOfMapChunksPerRow =
        (gridMapSize.mX + mapChunkWidth - 1) / mapChunkWidth;

    const int patchSize = gridMap.PatchSize();
    const int numOfPatchesX = gridMap.NumPatchX();
    const int numOfPatchesY = gridMap.NumPatchY();

    const Point2D<int> patchMinIdx = gridMap.CellIndexToPatchIndex(
        gridMapMinIdx.mX, gridMapMinIdx.mY);
    const Point2D<int> patchMinOffset = gridMap.CellIndexToPatchOffset(
        gridMapMinIdx.mX, gridMapMinIdx.mY);

    /* Initialize the index of the current patch being used */
    Point2D<int> patchIdx = patchMinIdx;
    Point2D<int> patchOffset = patchMinOffset;
    /* Initialize the pointer to the current patch */
    auto* pPatch = gridMap.PatchPtrAt(patchIdx);
    bool isPatchAllocated = pPatch->IsAllocated();

    for (int y = 0; y < gridMapSize.mY; ++y) {
        for (int chunkX = 0; chunkX < numOfMapChunksPerRow; ++chunkX) {
            /* Compute the grid map chunk (consecutive grid map elements) */
            std::uint64_t gridMapChunk = 0;

            for (int i = 0; i < mapChunkWidth; ++i) {
                /* Retrieve the occupancy probability value */
                const std::uint16_t rawValue = isPatchAllocated ?
                    pPatch->RawValue(patchOffset) : unknownRawValue;
                const std::uint64_t discretizedValue =
                    static_cast<std::uint64_t>(rawValue & 0xFF00);

                /* Pack the 8-bit discretized occupancy probability value to
                 * the grid map chunk (consecutive grid map elements) */
                gridMapChunk >>= 8;
                gridMapChunk |= (discretizedValue << 48);

                /* Update the index of the current patch */
                const bool isEndOfPatchX = patchOffset.mX == (patchSize - 1);
                patchOffset.mX = isEndOfPatchX ? 0 : patchOffset.mX + 1;
                patchIdx.mX = isEndOfPatchX ? patchIdx.mX + 1 : patchIdx.mX;
                patchIdx.mX = std::min(patchIdx.mX, numOfPatchesX - 1);

                /* Update the pointer to the current patch */
                pPatch = isEndOfPatchX ? gridMap.PatchPtrAt(patchIdx) : pPatch;
                /* Update the flag that indicates whether the current patch
                 * is allocated and has valid occupancy probability values */
                isPatchAllocated = isEndOfPatchX ?
                    pPatch->IsAllocated() : isPatchAllocated;
            }

            /* Write the grid map chunk */
            *pInput++ = gridMapChunk;
        }

        /* Update the index of the current patch */
        const bool isEndOfPatchY = patchOffset.mY == (patchSize - 1);
        patchOffset.mX = patchMinOffset.mX;
        patchOffset.mY = isEndOfPatchY ? 0 : patchOffset.mY + 1;
        patchIdx.mX = patchMinIdx.mX;
        patchIdx.mY = isEndOfPatchY ? patchIdx.mY + 1 : patchIdx.mY;
        patchIdx.mY = std::min(patchIdx.mY, numOfPatchesY - 1);

        /* Update the pointer to the current patch */
        pPatch = gridMap.PatchPtrAt(patchIdx);
        /* Update the flag that indicates whether the current patch
         * is allocated and has valid occupancy probability values */
        isPatchAllocated = pPatch->IsAllocated();
    }

    /* Transfer the grid map using the AXI DMA IP core */
    const std::size_t transferLengthInBytes =
        (gridMapSize.mY * numOfMapChunksPerRow + 1) * sizeof(std::uint64_t);
    this->mAxiDma[coreId]->SendChannel().Transfer(
        transferLengthInBytes, this->mInputData[coreId].PhysicalAddress());

    /* Wait for the data transfer to complete by polling the status register
     * of the AXI DMA IP core */
    if (this->mCommonConfig.mWaitForDmaTransfer)
        this->mAxiDma[coreId]->SendChannel().Wait();

    /* Update the metrics */
    const int numOfChunksTransferred = gridMapSize.mY * numOfMapChunksPerRow;
    this->mMetrics.mMapChunks->Increment(numOfChunksTransferred);
}

/* Receive the result through AXI DMA */
void ScanMatcherCorrelativeFPGA::ReceiveResult(
    const int coreId,
    int& scoreMax, int& bestX, int& bestY, int& bestTheta)
{
    /* Receive the result using the AXI DMA IP core */
    const std::size_t receiveLengthInBytes = 2 * sizeof(std::uint64_t);
    this->mAxiDma[coreId]->RecvChannel().Transfer(
        receiveLengthInBytes, this->mOutputData[coreId].PhysicalAddress());

    /* Wait for the data transfer to complete by polling the status register
     * of the AXI DMA IP core */
    if (this->mCommonConfig.mWaitForDmaTransfer)
        this->mAxiDma[coreId]->RecvChannel().Wait();

    /* Retrieve the pointer to the CMA memory */
    volatile std::uint64_t* pOutput =
        this->mOutputData->Ptr<volatile std::uint64_t>();

    /* Read the maximum score and the best X-coordinate */
    UnpackInt32(*pOutput++, scoreMax, bestX);
    UnpackInt32(*pOutput++, bestY, bestTheta);
}

/* Start the scan matcher IP core */
void ScanMatcherCorrelativeFPGA::StartIPCore(const int coreId)
{
    const std::uint32_t apStart = ToUnderlying(AxiLiteSApCtrl::Start);

    /* Set the control register of the scan matcher IP core */
    this->WriteCtrlReg(coreId, apStart);

    /* Wait until the scan matcher IP core starts by polling the
     * control registers of the AXI4-Lite slave interface */
    if (this->mCommonConfig.mWaitForCtrlReg)
        while (!(this->ReadCtrlReg(coreId) & apStart));
}

/* Wait until the scan matcher IP core is in idle state */
void ScanMatcherCorrelativeFPGA::WaitIPCore(const int coreId)
{
    const std::uint32_t apIdle = ToUnderlying(AxiLiteSApCtrl::Idle);

    /* Wait until the scan matcher IP core is in idle state by polling the
     * control registers of the AXI4-Lite slave interface */
    if (this->mCommonConfig.mWaitForCtrlReg)
        while (!(this->ReadCtrlReg(coreId) & apIdle));
}

} /* namespace MyGMapping */
} /* namespace Mapping */
