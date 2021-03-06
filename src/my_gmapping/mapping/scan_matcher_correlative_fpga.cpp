
/* scan_matcher_correlative_fpga.cpp */

#include "my_gmapping/mapping/scan_matcher_correlative_fpga.hpp"

#include <functional>
#include <future>

namespace MyGMapping {
namespace Mapping {

using namespace Hardware;

/* Static assertions for safety */
static_assert(std::is_same<GridMap::GridType::ValueType, std::uint16_t>::value,
              "Grid map should store occupancy probability values as "
              "discretized 16-bit unsigned integers");
static_assert(GridMap::GridType::UnknownValue == 0,
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
    mMapTransfers(nullptr),
    mMapTransferSkips(nullptr),
    mMapChunks(nullptr),
    mScanTransfers(nullptr),
    mScanTransferSkips(nullptr),
    mNumOfTransferredScans(nullptr),
    mScoreValue(nullptr),
    mLikelihoodValue(nullptr)
{
    /* Retrieve the metrics manager instance */
    auto* const pMetricManager = Metric::MetricManager::Instance();

    /* Register the value sequence metrics */
    this->mInputSetupTime = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".InputSetupTime");
    this->mSetupIPTime = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".SetupIPTime");
    this->mScanSendTime = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".ScanSendTime");
    this->mMapSendTime = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".MapSendTime");
    this->mOptimizationTime = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".OptimizationTime");
    this->mWaitIPTime = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".WaitIPTime");
    this->mScanMatchingTime = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".ScanMatchingTime");

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

    this->mMapTransfers = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".MapTransfers");
    this->mMapTransferSkips = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".MapTransferSkips");
    this->mMapChunks = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".MapChunks");
    this->mScanTransfers = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".ScanTransfers");
    this->mScanTransferSkips = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".ScanTransferSkips");
    this->mNumOfTransferredScans = pMetricManager->AddValueSequence<int>(
        scanMatcherName + ".NumOfTransferredScans");

    this->mScoreValue = pMetricManager->AddValueSequence<float>(
        scanMatcherName + ".ScoreValue");
    this->mLikelihoodValue = pMetricManager->AddValueSequence<float>(
        scanMatcherName + ".LikelihoodValue");
}

/* Reserve the buffer to store the processing times */
void ScanMatcherCorrelativeFPGAMetrics::Resize(
    const std::size_t numOfParticles)
{
    this->mTimes.resize(numOfParticles);
    this->mParams.resize(numOfParticles);
}

/* Set the processing times for each particle */
void ScanMatcherCorrelativeFPGAMetrics::SetTimes(
    const std::size_t idx, const Times& times)
{
    this->mTimes[idx] = times;
}

/* Set the parameter settings for each particle */
void ScanMatcherCorrelativeFPGAMetrics::SetParameters(
    const std::size_t idx, const Parameters& parameters)
{
    this->mParams[idx] = parameters;
}

/* Collect the particle-wise metrics and update the overall metrics */
void ScanMatcherCorrelativeFPGAMetrics::Update(
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
    this->mSetupIPTime->Observe(collectTimes(
        [](int value, const Times& times) {
            return value + times.mSetupIPTime; }));
    this->mScanSendTime->Observe(collectTimes(
        [](int value, const Times& times) {
            return value + times.mScanSendTime; }));
    this->mMapSendTime->Observe(collectTimes(
        [](int value, const Times& times) {
            return value + times.mMapSendTime; }));
    this->mOptimizationTime->Observe(collectTimes(
        [](int value, const Times& times) {
            return value + times.mOptimizationTime; }));
    this->mWaitIPTime->Observe(collectTimes(
        [](int value, const Times& times) {
            return value + times.mWaitIPTime; }));
    this->mScanMatchingTime->Observe(collectTimes(
        [](int value, const Times& times) {
            return value + times.mScanMatchingTime; }));

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

    /* Store the collected metric values */
    this->mMapTransfers->Observe(collectParams(
        [](int value, const Parameters& params) {
            return params.mMapTransferred ? (value + 1) : value; }));
    this->mMapTransferSkips->Observe(collectParams(
        [](int value, const Parameters& params) {
            return params.mMapTransferred ? value : (value + 1); }));
    this->mMapChunks->Observe(collectParams(
        [](int value, const Parameters& params) {
            return value + params.mMapChunks; }));
    this->mScanTransfers->Observe(collectParams(
        [](int value, const Parameters& params) {
            return params.mScanTransferred ? (value + 1) : value; }));
    this->mScanTransferSkips->Observe(collectParams(
        [](int value, const Parameters& params) {
            return params.mScanTransferred ? value : (value + 1); }));
    this->mNumOfTransferredScans->Observe(collectParams(
        [](int value, const Parameters& params) {
            return value + params.mNumOfTransferredScans; }));
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
    XAssert(this->mCommonConfig.mMapChunkWidth == 8,
            "Width of the grid map chunk (consecutive grid map cells) "
            "should be exactly 8");

    /* Initialize the scan matcher IP cores */
    this->Initialize(0);
    this->Initialize(1);
}

/* Optimize the particle poses based on the correlative scan matching */
ScanMatchingResultVector ScanMatcherCorrelativeFPGA::OptimizePose(
    const ScanMatchingQueryVector& queries,
    const Sensor::ScanDataPtr<double>& scanData)
{
    /* Check the number of the particles */
    XAssert(queries.size() % NumOfIPCores == 0,
            "Number of the particles must be multiple of the "
            "number of the scan matcher IP cores implemented on the device");

    /* Resize the buffer to store the metrics for each particle */
    this->mMetrics.Resize(queries.size());

    const std::size_t size = queries.size();
    const std::size_t half = queries.size() / 2;

    /* Start a sub-thread to process the second half of the particles */
    auto futSubResults = std::async(std::launch::async, [&]() {
        return this->OptimizePoseCore(1, half, size, queries, scanData); });
    /* Process the first half of the particles */
    auto results = this->OptimizePoseCore(0, 0, half, queries, scanData);
    /* Retrieve the second half of the results */
    auto subResults = futSubResults.get();

    /* Concatenate the results */
    results.insert(results.end(),
                   std::make_move_iterator(subResults.begin()),
                   std::make_move_iterator(subResults.end()));

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

/* Optimize the particle poses using the scan matcher IP core */
ScanMatchingResultVector ScanMatcherCorrelativeFPGA::OptimizePoseCore(
    const int coreId,
    const std::size_t idxBegin,
    const std::size_t idxEnd,
    const ScanMatchingQueryVector& queries,
    const Sensor::ScanDataPtr<double>& scanData)
{
    ScanMatchingResultVector results;
    results.reserve(idxEnd - idxBegin);

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
    for (std::size_t i = idxBegin; i < idxEnd; ++i) {
        /* Create the timer */
        Metric::Timer outerTimer;
        Metric::Timer timer;
        /* Metric values for each particle */
        ScanMatcherMetrics::Times times;
        ScanMatcherMetrics::Parameters params;

        const auto& initialPose = queries[i].mInitialPose;
        const auto& gridMap = queries[i].mGridMap;

        /* Compute the sensor pose from the initial particle pose */
        const RobotPose2D<double> sensorPose =
            Compound(initialPose, scanData->RelativeSensorPose());
        /* Compute the minimum possible sensor pose */
        const RobotPose2D<double> minSensorPose {
            sensorPose.mX - stepX * winX,
            sensorPose.mY - stepY * winY,
            sensorPose.mTheta - stepTheta * winTheta };

        /* Compute the bounding box of the grid map */
        const BoundingBox<int> boundingBox =
            this->ComputeBoundingBox(gridMap, sensorPose);
        /* Compute the size of the cropped grid map */
        const Point2D<int> gridMapSize {
            boundingBox.Width(), boundingBox.Height() };
        /* Compute the minimum position of the cropped grid map */
        const Point2D<double> gridMapMinPos =
            gridMap.IndexToPosition(boundingBox.mMin.mY, boundingBox.mMin.mX);

        /* Update the metrics and restart the timer */
        times.mInputSetupTime = timer.ElapsedMicro();
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
        times.mSetupIPTime = timer.ElapsedMicro();
        timer.Start();

        /* Send the scan data for the first particle only */
        const bool scanDataTransferred = (i == idxBegin);
        this->SendScanData(coreId, scanDataTransferred, scanData, params);
        /* Update the metrics and restart the timer */
        times.mScanSendTime = timer.ElapsedMicro();
        timer.Start();

        /* Send the grid map */
        this->SendGridMap(coreId, gridMap, boundingBox, params);
        /* Update the metrics and restart the timer */
        times.mMapSendTime = timer.ElapsedMicro();
        timer.Start();

        /* Receive the result */
        int scoreMax;
        int bestX;
        int bestY;
        int bestTheta;
        this->ReceiveResult(coreId, scoreMax, bestX, bestY, bestTheta);
        /* Update the metrics and restart the timer */
        times.mOptimizationTime = timer.ElapsedMicro();
        timer.Start();

        /* Wait for the scan matcher IP core */
        this->WaitIPCore(coreId);
        /* Update the metrics and stop the timer */
        times.mWaitIPTime = timer.ElapsedMicro();
        timer.Stop();

        /* Compute the best sensor pose */
        const RobotPose2D<double> bestSensorPose {
            minSensorPose.mX + bestX * stepX,
            minSensorPose.mY + bestY * stepY,
            minSensorPose.mTheta + bestTheta * stepTheta };
        /* Compute the estimated robot pose and the likelihood value */
        const RobotPose2D<double> estimatedPose =
            MoveBackward(bestSensorPose, scanData->RelativeSensorPose());
        /* Compute the difference from the initial pose */
        const RobotPose2D<double> diffPose =
            InverseCompound(initialPose, estimatedPose);

        /* Set the resulting score */
        const double factor = (1 << this->mCommonConfig.mMapBitWidth) - 1;
        const double score = static_cast<double>(scoreMax) / factor;
        const double normalizedScore = score / numOfScansTransferred;

        /* Set the resulting likelihood */
        const double likelihood = this->mLikelihoodFunc->Likelihood(
            gridMap, scanData, bestSensorPose);
        const double normalizedLikelihood = likelihood / scanData->NumOfScans();

        /* Update the metrics */
        times.mScanMatchingTime = outerTimer.ElapsedMicro();
        outerTimer.Stop();

        params.mDiffTranslation = Distance(diffPose);
        params.mDiffRotation = std::abs(diffPose.mTheta);
        params.mWinSizeX = winX;
        params.mWinSizeY = winY;
        params.mWinSizeTheta = winTheta;
        params.mStepSizeX = stepX;
        params.mStepSizeY = stepY;
        params.mStepSizeTheta = stepTheta;
        params.mScoreValue = normalizedScore;
        params.mLikelihoodValue = normalizedLikelihood;

        /* Set the metrics */
        this->mMetrics.SetTimes(i, times);
        this->mMetrics.SetParameters(i, params);

        /* Append the scan matching result */
        results.emplace_back(initialPose, estimatedPose,
                             normalizedLikelihood, likelihood,
                             normalizedScore, score);
    }

    return results;
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

/* Compute the bounding box of the grid map */
BoundingBox<int> ScanMatcherCorrelativeFPGA::ComputeBoundingBox(
    const GridMap& gridMap,
    const RobotPose2D<double>& sensorPose) const
{
    /* Make sure that the block size is the multiple of 4 so that the
     * size of the grid map `gridMap.Cols()` and `gridMap.Rows()` are
     * also the multiples of 4 */
    Assert(gridMap.BlockSize() % 4 == 0);

    const int mapCols = gridMap.Cols();
    const int mapRows = gridMap.Rows();

    /* Desired size of the grid map */
    const int mapColsMax = this->mCommonConfig.mMaxMapSizeX;
    const int mapRowsMax = this->mCommonConfig.mMaxMapSizeY;

    /* Compute the center position of the cropped grid map */
    const Point2D<int> desiredCenterIdx =
        gridMap.PositionToIndex(sensorPose.mX, sensorPose.mY);
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
    const auto alignBy4 = [](const int x) { return (x >> 2) << 2; };
    const int colMin = alignBy4(centerIdx.mX - mapColsMax / 2);
    const int rowMin = alignBy4(centerIdx.mY - mapRowsMax / 2);
    const int colMax = alignBy4(centerIdx.mX + mapColsMax / 2);
    const int rowMax = alignBy4(centerIdx.mY + mapRowsMax / 2);

    const Point2D<int> idxMin { std::max(colMin, 0),
                                std::max(rowMin, 0) };
    const Point2D<int> idxMax { std::min(colMax, gridMap.Cols()),
                                std::min(rowMax, gridMap.Rows()) };

    /* Return the bounding box */
    return BoundingBox<int> { idxMin, idxMax };
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
    const auto doubleToU32 = [](const double x) {
        return FloatToU32(static_cast<float>(x)); };
    const auto writeReg = [&](const std::uint32_t offsetInBytes,
                              const std::uint32_t value) {
        this->mControlRegisters[coreId]->Write(offsetInBytes, value); };
    const auto& config = this->mCommonConfig;

    /* Set the actual number of the scan points */
    writeReg(config.mAxiLiteSNumOfScans, I32ToU32(numOfScans));
    /* Set the maximum scan range considered valid */
    writeReg(config.mAxiLiteSScanRangeMax, doubleToU32(scanRangeMax));
    /* Set the score threshold (for loop detection) */
    writeReg(config.mAxiLiteSScoreThreshold, I32ToU32(scoreThreshold));

    /* Set the minimum possible sensor pose */
    writeReg(config.mAxiLiteSPoseX, doubleToU32(minSensorPose.mX));
    writeReg(config.mAxiLiteSPoseY, doubleToU32(minSensorPose.mY));
    writeReg(config.mAxiLiteSPoseTheta, doubleToU32(minSensorPose.mTheta));

    /* Set the actual size of the cropped grid map */
    writeReg(config.mAxiLiteSMapSizeX, I32ToU32(gridMapSize.mX));
    writeReg(config.mAxiLiteSMapSizeY, I32ToU32(gridMapSize.mY));
    /* Set the minimum coordinate of the cropped grid map */
    writeReg(config.mAxiLiteSMapMinX, doubleToU32(gridMapMinPos.mX));
    writeReg(config.mAxiLiteSMapMinY, doubleToU32(gridMapMinPos.mY));

    /* Set the size of the search window */
    writeReg(config.mAxiLiteSWinX, I32ToU32(winX));
    writeReg(config.mAxiLiteSWinY, I32ToU32(winY));
    writeReg(config.mAxiLiteSWinTheta, I32ToU32(winTheta));
    /* Set the search step */
    writeReg(config.mAxiLiteSStepX, doubleToU32(stepX));
    writeReg(config.mAxiLiteSStepY, doubleToU32(stepY));
    writeReg(config.mAxiLiteSStepTheta, doubleToU32(stepTheta));
}

/* Send the scan data through AXI DMA */
void ScanMatcherCorrelativeFPGA::SendScanData(
    const int coreId,
    const bool scanDataTransferred,
    const Sensor::ScanDataPtr<double>& scanData,
    ScanMatcherMetrics::Parameters& params)
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
        this->mAxiDma[coreId]->SendChannel().Wait();

        /* Update the metrics */
        params.mScanTransferred = false;
        params.mNumOfTransferredScans = 0;

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
    this->mAxiDma[coreId]->SendChannel().Wait();

    /* Update the metrics */
    params.mScanTransferred = true;
    params.mNumOfTransferredScans = numOfScansTransferred;
}

/* Send the grid map through AXI DMA */
void ScanMatcherCorrelativeFPGA::SendGridMap(
    const int coreId,
    const GridMap& gridMap,
    const BoundingBox<int>& desiredBox,
    ScanMatcherMetrics::Parameters& params)
{
    /* Retrieve the pointer to the CMA memory */
    volatile std::uint64_t* pInput =
        this->mInputData[coreId].Ptr<volatile std::uint64_t>();

    /* Write the flag to indicate that the grid map is transferred */
    *pInput++ = static_cast<std::uint64_t>(1);

    /* Write the cropped grid map */
    const int chunkWidth = this->mCommonConfig.mMapChunkWidth;
    const int chunkCols = (desiredBox.Width() + chunkWidth - 1) / chunkWidth;
    const BoundingBox<int> boundingBox {
        desiredBox.mMin.mX, desiredBox.mMin.mY,
        desiredBox.mMin.mX + chunkCols * chunkWidth,
        desiredBox.mMin.mY + desiredBox.Height() };

    auto* pBuffer = this->mInputData[coreId].Ptr<std::uint32_t>() +
                    sizeof(std::uint64_t) / sizeof(std::uint32_t);
    gridMap.CopyValuesU8x4(pBuffer, boundingBox);

    /* Transfer the grid map using the AXI DMA IP core */
    const std::size_t transferLengthInBytes =
        (desiredBox.Height() * chunkCols + 1) * sizeof(std::uint64_t);
    this->mAxiDma[coreId]->SendChannel().Transfer(
        transferLengthInBytes, this->mInputData[coreId].PhysicalAddress());
    this->mAxiDma[coreId]->SendChannel().Wait();

    /* Update the metrics */
    params.mMapTransferred = true;
    params.mMapChunks = desiredBox.Height() * chunkCols;
}

/* Receive the result through AXI DMA */
void ScanMatcherCorrelativeFPGA::ReceiveResult(
    const int coreId, int& scoreMax, int& bestX, int& bestY, int& bestTheta)
{
    /* Receive the result using the AXI DMA IP core */
    const std::size_t receiveLengthInBytes = 2 * sizeof(std::uint64_t);
    this->mAxiDma[coreId]->RecvChannel().Transfer(
        receiveLengthInBytes, this->mOutputData[coreId].PhysicalAddress());
    this->mAxiDma[coreId]->RecvChannel().Wait();

    /* Retrieve the pointer to the CMA memory */
    volatile std::uint64_t* pOutput =
        this->mOutputData->Ptr<volatile std::uint64_t>();

    /* Read the maximum score and the best X-coordinate */
    UnpackI32(*pOutput++, scoreMax, bestX);
    UnpackI32(*pOutput++, bestY, bestTheta);
}

/* Start the scan matcher IP core */
void ScanMatcherCorrelativeFPGA::StartIPCore(const int coreId)
{
    /* Set the control register of the scan matcher IP core */
    this->WriteCtrlReg(coreId, ToUnderlying(AxiLiteSApCtrl::Start));
    /* Wait until the scan matcher IP core starts by polling the
     * control registers of the AXI4-Lite slave interface */
    while (!(this->ReadCtrlReg(coreId) & ToUnderlying(AxiLiteSApCtrl::Start)));
}

/* Wait until the scan matcher IP core is in idle state */
void ScanMatcherCorrelativeFPGA::WaitIPCore(const int coreId)
{
    /* Wait until the scan matcher IP core is in idle state by polling the
     * control registers of the AXI4-Lite slave interface */
    while (!(this->ReadCtrlReg(coreId) & ToUnderlying(AxiLiteSApCtrl::Idle)));
}

} /* namespace MyGMapping */
} /* namespace Mapping */
