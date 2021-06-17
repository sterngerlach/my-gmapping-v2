
/* slam_launcher.cpp */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>

#ifdef __GNUC__
#if (__GNUC__ >= 6) && (__GNUC__ < 8)
#include <experimental/filesystem>
#elif (__GNUC__ >= 8)
#include <filesystem>
#endif
#endif

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "my_gmapping/hw/bitstream_loader.hpp"
#include "my_gmapping/hw/cma_manager.hpp"
#include "my_gmapping/io/gnuplot_helper.hpp"
#include "my_gmapping/io/map_saver.hpp"
#include "my_gmapping/io/carmen/carmen_reader.hpp"
#include "my_gmapping/mapping/covariance_estimator.hpp"
#include "my_gmapping/mapping/grid_map_builder.hpp"
#include "my_gmapping/mapping/likelihood_greedy_endpoint.hpp"
#include "my_gmapping/mapping/map_builder.hpp"
#include "my_gmapping/mapping/motion_model_relative_pose.hpp"
#include "my_gmapping/mapping/scan_interpolator.hpp"
#include "my_gmapping/mapping/scan_matcher.hpp"
#include "my_gmapping/mapping/scan_matcher_correlative.hpp"
#include "my_gmapping/mapping/scan_matcher_correlative_fpga.hpp"
#include "my_gmapping/mapping/scan_matcher_hill_climbing.hpp"
#include "my_gmapping/mapping/weight_normalizer.hpp"
#include "my_gmapping/metric/metric.hpp"
#include "my_gmapping/sensor/sensor_data.hpp"

using namespace MyGMapping;

/* Declare namespaces for convenience */
namespace pt = boost::property_tree;

#ifdef __GNUC__
#if (__GNUC__ >= 6) && (__GNUC__ < 8)
namespace fs = std::experimental::filesystem;
#elif (__GNUC__ >= 8)
namespace fs = std::filesystem;
#endif
#endif

/* Create the robot motion model */
std::unique_ptr<Mapping::MotionModel> CreateMotionModel(
    const pt::ptree& jsonSettings)
{
    /* Load settings for the robot motion model */
    const double sigmaLinear = jsonSettings.get<double>(
        "MotionModelRelativePose.SigmaLinear");
    const double sigmaAngular = jsonSettings.get<double>(
        "MotionModelRelativePose.SigmaAngular");
    const double sigmaLinearToAngular = jsonSettings.get<double>(
        "MotionModelRelativePose.SigmaLinearToAngular");
    const double sigmaAngularToLinear = jsonSettings.get<double>(
        "MotionModelRelativePose.SigmaAngularToLinear");

    /* Create the robot motion model */
    auto pMotionModel = std::make_unique<Mapping::MotionModelRelativePose>(
        sigmaLinear, sigmaAngular,
        sigmaLinearToAngular, sigmaAngularToLinear);

    return pMotionModel;
}

/* Create the greedy endpoint likelihood function */
std::unique_ptr<Mapping::LikelihoodFunction> CreateLikelihoodGreedyEndpoint(
    const pt::ptree& jsonSettings)
{
    /* Load settings for the greedy endpoint likelihood function */
    const double mapResolution = jsonSettings.get<double>(
        "LikelihoodGreedyEndpoint.MapResolution");
    const double minUsableRange = jsonSettings.get<double>(
        "LikelihoodGreedyEndpoint.MinUsableRange");
    const double maxUsableRange = jsonSettings.get<double>(
        "LikelihoodGreedyEndpoint.MaxUsableRange");
    const double hitAndMissedCellDist = jsonSettings.get<double>(
        "LikelihoodGreedyEndpoint.HitAndMissedCellDist");
    const double occupancyThreshold = jsonSettings.get<double>(
        "LikelihoodGreedyEndpoint.OccupancyThreshold");
    const double gaussianSigma = jsonSettings.get<double>(
        "LikelihoodGreedyEndpoint.GaussianSigma");
    const int kernelSize = jsonSettings.get<int>(
        "LikelihoodGreedyEndpoint.KernelSize");
    const double likelihoodScale = jsonSettings.get<double>(
        "LikelihoodGreedyEndpoint.LikelihoodScale");

    /* Create the greedy endpoint likelihood function */
    auto pLikelihoodFunc = std::make_unique<Mapping::LikelihoodGreedyEndpoint>(
        mapResolution, minUsableRange, maxUsableRange, hitAndMissedCellDist,
        occupancyThreshold, gaussianSigma, kernelSize, likelihoodScale);

    return pLikelihoodFunc;
}

/* Create the hill-climbing based scan matcher */
std::unique_ptr<Mapping::ScanMatcher> CreateScanMatcherHillClimbing(
    const pt::ptree& jsonSettings,
    std::unique_ptr<Mapping::LikelihoodFunction>&& pLikelihoodFunc)
{
    /* Load settings for the hill-climbing based scan matcher */
    const std::string scanMatcherName = jsonSettings.get<std::string>(
        "ScanMatcherHillClimbing.Name");
    const double linearDelta = jsonSettings.get<double>(
        "ScanMatcherHillClimbing.LinearDelta");
    const double angularDelta = jsonSettings.get<double>(
        "ScanMatcherHillClimbing.AngularDelta");
    const int maxIterations = jsonSettings.get<int>(
        "ScanMatcherHillClimbing.MaxIterations");
    const int numOfRefinements = jsonSettings.get<int>(
        "ScanMatcherHillClimbing.NumOfRefinements");

    /* Create the hill-climbing based scan matcher */
    auto pScanMatcher = std::make_unique<Mapping::ScanMatcherHillClimbing>(
        scanMatcherName, std::move(pLikelihoodFunc),
        linearDelta, angularDelta, maxIterations, numOfRefinements);

    return pScanMatcher;
}

/* Create the real-time correlative-based scan matcher */
std::unique_ptr<Mapping::ScanMatcher> CreateScanMatcherCorrelative(
    const pt::ptree& jsonSettings,
    std::unique_ptr<Mapping::LikelihoodFunction>&& pLikelihoodFunc)
{
    /* Load settings for the real-time correlative-based scan matcher */
    const std::string scanMatcherName = jsonSettings.get<std::string>(
        "ScanMatcherCorrelative.Name");
    const bool useCroppedMap = jsonSettings.get<bool>(
        "ScanMatcherCorrelative.UseCroppedMap");
    const int croppedMapSizeX = jsonSettings.get<int>(
        "ScanMatcherCorrelative.CroppedMapSizeX");
    const int croppedMapSizeY = jsonSettings.get<int>(
        "ScanMatcherCorrelative.CroppedMapSizeY");
    const int lowResolution = jsonSettings.get<int>(
        "ScanMatcherCorrelative.LowResolution");
    const double rangeX = jsonSettings.get<double>(
        "ScanMatcherCorrelative.RangeX");
    const double rangeY = jsonSettings.get<double>(
        "ScanMatcherCorrelative.RangeY");
    const double rangeTheta = jsonSettings.get<double>(
        "ScanMatcherCorrelative.RangeTheta");
    const double scanRangeMax = jsonSettings.get<double>(
        "ScanMatcherCorrelative.ScanRangeMax");

    /* Create the real-time correlative-based scan matcher */
    auto pScanMatcher = std::make_unique<Mapping::ScanMatcherCorrelative>(
        scanMatcherName, std::move(pLikelihoodFunc),
        useCroppedMap, croppedMapSizeX, croppedMapSizeY,
        lowResolution, rangeX, rangeY, rangeTheta, scanRangeMax);

    return pScanMatcher;
}

/* Create the real-time correlative-based scan matcher IP core interface */
std::unique_ptr<Mapping::ScanMatcher> CreateScanMatcherCorrelativeFPGA(
    const pt::ptree& jsonSettings,
    std::unique_ptr<Mapping::LikelihoodFunction>&& pLikelihoodFunc)
{
    /* Load settings for the real-time correlative-based scan matcher */
    Mapping::ScanMatcherCorrelativeFPGACommonConfig commonConfig;
    Mapping::ScanMatcherCorrelativeFPGAConfig scanMatcherConfig0;
    Mapping::ScanMatcherCorrelativeFPGAConfig scanMatcherConfig1;
    Mapping::AxiDmaConfig axiDmaConfig0;
    Mapping::AxiDmaConfig axiDmaConfig1;

    const pt::ptree& scanMatcherSettings =
        jsonSettings.get_child("ScanMatcherCorrelativeFPGA");

    /* Load settings for the register offsets */
    const pt::ptree& registerOffsetSettings =
        scanMatcherSettings.get_child("RegisterOffsets");

    const auto offsetString = [&](const std::string& registerName) {
        return registerOffsetSettings.get<std::string>(registerName); };

    const auto offsetApCtrl = offsetString("AxiLiteSApCtrl");
    const auto offsetGIE = offsetString("AxiLiteSGIE");
    const auto offsetIER = offsetString("AxiLiteSIER");
    const auto offsetISR = offsetString("AxiLiteSISR");
    const auto offsetNumOfScans = offsetString("AxiLiteSNumOfScans");
    const auto offsetScanRangeMax = offsetString("AxiLiteSScanRangeMax");
    const auto offsetScoreThreshold = offsetString("AxiLiteSScoreThreshold");
    const auto offsetPoseX = offsetString("AxiLiteSPoseX");
    const auto offsetPoseY = offsetString("AxiLiteSPoseY");
    const auto offsetPoseTheta = offsetString("AxiLiteSPoseTheta");
    const auto offsetMapSizeX = offsetString("AxiLiteSMapSizeX");
    const auto offsetMapSizeY = offsetString("AxiLiteSMapSizeY");
    const auto offsetMapMinX = offsetString("AxiLiteSMapMinX");
    const auto offsetMapMinY = offsetString("AxiLiteSMapMinY");
    const auto offsetWinX = offsetString("AxiLiteSWinX");
    const auto offsetWinY = offsetString("AxiLiteSWinY");
    const auto offsetWinTheta = offsetString("AxiLiteSWinTheta");
    const auto offsetStepX = offsetString("AxiLiteSStepX");
    const auto offsetStepY = offsetString("AxiLiteSStepY");
    const auto offsetStepTheta = offsetString("AxiLiteSStepTheta");

    /* Set the register offsets */
    const auto toAddress = [](const std::string& value) {
        return std::stoul(value, nullptr, 0); };

    commonConfig.mAxiLiteSApCtrl = toAddress(offsetApCtrl);
    commonConfig.mAxiLiteSGIE = toAddress(offsetGIE);
    commonConfig.mAxiLiteSIER = toAddress(offsetIER);
    commonConfig.mAxiLiteSISR = toAddress(offsetISR);
    commonConfig.mAxiLiteSNumOfScans = toAddress(offsetNumOfScans);
    commonConfig.mAxiLiteSScanRangeMax = toAddress(offsetScanRangeMax);
    commonConfig.mAxiLiteSScoreThreshold = toAddress(offsetScoreThreshold);
    commonConfig.mAxiLiteSPoseX = toAddress(offsetPoseX);
    commonConfig.mAxiLiteSPoseY = toAddress(offsetPoseY);
    commonConfig.mAxiLiteSPoseTheta = toAddress(offsetPoseTheta);
    commonConfig.mAxiLiteSMapSizeX = toAddress(offsetMapSizeX);
    commonConfig.mAxiLiteSMapSizeY = toAddress(offsetMapSizeY);
    commonConfig.mAxiLiteSMapMinX = toAddress(offsetMapMinX);
    commonConfig.mAxiLiteSMapMinY = toAddress(offsetMapMinY);
    commonConfig.mAxiLiteSWinX = toAddress(offsetWinX);
    commonConfig.mAxiLiteSWinY = toAddress(offsetWinY);
    commonConfig.mAxiLiteSWinTheta = toAddress(offsetWinTheta);
    commonConfig.mAxiLiteSStepX = toAddress(offsetStepX);
    commonConfig.mAxiLiteSStepY = toAddress(offsetStepY);
    commonConfig.mAxiLiteSStepTheta = toAddress(offsetStepTheta);

    /* Load the hardware-specific parameters */
    const auto getString = [&](const std::string& configName) {
        return scanMatcherSettings.get<std::string>(configName); };
    const auto getInt = [&](const std::string& configName) {
        return scanMatcherSettings.get<int>(configName); };
    const auto getDouble = [&](const std::string& configName) {
        return scanMatcherSettings.get<double>(configName); };

    const auto scanMatcherName = getString("Name");
    const auto maxNumOfScans = getInt("MaxNumOfScans");
    const auto mapResolution = getDouble("MapResolution");
    const auto maxMapSizeX = getInt("MaxMapSizeX");
    const auto maxMapSizeY = getInt("MaxMapSizeY");
    const auto lowResolution = getInt("CoarseMapResolution");
    const auto mapBitWidth = getInt("MapBitWidth");
    const auto mapChunkWidth = getInt("MapChunkWidth");

    commonConfig.mMaxNumOfScans = maxNumOfScans;
    commonConfig.mMapResolution = mapResolution;
    commonConfig.mMaxMapSizeX = maxMapSizeX;
    commonConfig.mMaxMapSizeY = maxMapSizeY;
    commonConfig.mLowResolution = lowResolution;
    commonConfig.mMapBitWidth = mapBitWidth;
    commonConfig.mMapChunkWidth = mapChunkWidth;

    const auto axiLiteSBaseAddress0 = getString("AxiLiteSBaseAddress0");
    const auto axiLiteSAddressRange0 = getString("AxiLiteSAddressRange0");
    const auto axiDmaBaseAddress0 = getString("AxiDmaBaseAddress0");
    const auto axiDmaAddressRange0 = getString("AxiDmaAddressRange0");

    const auto axiLiteSBaseAddress1 = getString("AxiLiteSBaseAddress1");
    const auto axiLiteSAddressRange1 = getString("AxiLiteSAddressRange1");
    const auto axiDmaBaseAddress1 = getString("AxiDmaBaseAddress1");
    const auto axiDmaAddressRange1 = getString("AxiDmaAddressRange1");

    scanMatcherConfig0.mAxiLiteSBaseAddress = toAddress(axiLiteSBaseAddress0);
    scanMatcherConfig0.mAxiLiteSAddressRange = toAddress(axiLiteSAddressRange0);
    axiDmaConfig0.mBaseAddress = toAddress(axiDmaBaseAddress0);
    axiDmaConfig0.mAddressRange = toAddress(axiDmaAddressRange0);

    scanMatcherConfig1.mAxiLiteSBaseAddress = toAddress(axiLiteSBaseAddress1);
    scanMatcherConfig1.mAxiLiteSAddressRange = toAddress(axiLiteSAddressRange1);
    axiDmaConfig1.mBaseAddress = toAddress(axiDmaBaseAddress1);
    axiDmaConfig1.mAddressRange = toAddress(axiDmaAddressRange1);

    /* Load the scan matching parameters */
    const auto searchRangeX = getDouble("SearchRangeX");
    const auto searchRangeY = getDouble("SearchRangeY");
    const auto searchRangeTheta = getDouble("SearchRangeTheta");
    const auto scanRangeMax = getDouble("ScanRangeMax");

    /* Create the real-time correlative-based scan matcher IP core interface */
    auto pScanMatcher = std::make_unique<Mapping::ScanMatcherCorrelativeFPGA>(
        scanMatcherName, std::move(pLikelihoodFunc),
        searchRangeX, searchRangeY, searchRangeTheta, scanRangeMax,
        std::move(scanMatcherConfig0), std::move(scanMatcherConfig1),
        std::move(commonConfig),
        std::move(axiDmaConfig0), std::move(axiDmaConfig1));

    return pScanMatcher;
}

/* Create the software scan matcher */
std::unique_ptr<Mapping::ScanMatcher> CreateScanMatcher(
    const pt::ptree& jsonSettings,
    const std::string& scanMatcherType,
    std::unique_ptr<Mapping::LikelihoodFunction>&& pLikelihoodFunc)
{
    if (scanMatcherType == "HillClimbing")
        return CreateScanMatcherHillClimbing(
            jsonSettings, std::move(pLikelihoodFunc));
    else if (scanMatcherType == "RealTimeCorrelative")
        return CreateScanMatcherCorrelative(
            jsonSettings, std::move(pLikelihoodFunc));

    return nullptr;
}

/* Create the hill-climbing based scan matcher */
std::unique_ptr<Mapping::ScanMatcher> CreateFinalScanMatcherHillClimbing(
    const pt::ptree& jsonSettings,
    std::unique_ptr<Mapping::LikelihoodFunction>&& pLikelihoodFunc)
{
    /* Load settings for the hill-climbing based scan matcher */
    const std::string scanMatcherName = jsonSettings.get<std::string>(
        "FinalScanMatcherHillClimbing.Name");
    const double linearDelta = jsonSettings.get<double>(
        "FinalScanMatcherHillClimbing.LinearDelta");
    const double angularDelta = jsonSettings.get<double>(
        "FinalScanMatcherHillClimbing.AngularDelta");
    const int maxIterations = jsonSettings.get<int>(
        "FinalScanMatcherHillClimbing.MaxIterations");
    const int numOfRefinements = jsonSettings.get<int>(
        "FinalScanMatcherHillClimbing.NumOfRefinements");

    /* Create the hill-climbing based scan matcher */
    auto pScanMatcher = std::make_unique<Mapping::ScanMatcherHillClimbing>(
        scanMatcherName, std::move(pLikelihoodFunc),
        linearDelta, angularDelta, maxIterations, numOfRefinements);

    return pScanMatcher;
}

/* Create the scan interpolator */
std::unique_ptr<Mapping::ScanInterpolator> CreateScanInterpolator(
    const pt::ptree& jsonSettings)
{
    /* Load settings for the scan interpolator */
    const double distScans = jsonSettings.get<double>(
        "ScanInterpolator.DistScans");
    const double distThresholdEmpty = jsonSettings.get<double>(
        "ScanInterpolator.DistThresholdEmpty");

    /* Create the scan interpolator */
    auto pScanInterpolator = std::make_unique<Mapping::ScanInterpolator>(
        distScans, distThresholdEmpty);

    return pScanInterpolator;
}

/* Create the pose covariance estimator */
std::unique_ptr<Mapping::CovarianceEstimator> CreateCovarianceEstimator(
    const pt::ptree& jsonSettings)
{
    /* Load settings for the pose covariance estimator */
    const double covarianceScale = jsonSettings.get<double>(
        "CovarianceEstimator.CovarianceScale");

    /* Create the pose covariance estimator */
    auto pCovarianceEstimator =
        std::make_unique<Mapping::CovarianceEstimator>(covarianceScale);

    return pCovarianceEstimator;
}

/* Determine weight normalization method */
Mapping::WeightNormalizationType DetermineWeightNormalizationType(
    const std::string& weightType)
{
    if (weightType == "Basic")
        return Mapping::WeightNormalizationType::Basic;
    else if (weightType == "ExponentialWeight")
        return Mapping::WeightNormalizationType::ExponentialWeight;

    /* Return the default normalization method */
    return Mapping::WeightNormalizationType::Basic;
}

/* Create the grid map builder */
std::shared_ptr<Mapping::GridMapBuilder> CreateGridMapBuilder(
    const pt::ptree& jsonSettings)
{
    /* Create the robot motion model */
    auto pMotionModel = CreateMotionModel(jsonSettings);

    /* Create the likelihood function for the grid map builder */
    auto pLikelihoodFunc = CreateLikelihoodGreedyEndpoint(jsonSettings);

    /* Load settings for hardware offloading */
    const bool useHardwareScanMatcher = jsonSettings.get<bool>(
        "GridMapBuilder.UseHardwareScanMatcher");

    std::unique_ptr<Mapping::ScanMatcher> pScanMatcher;

    if (useHardwareScanMatcher) {
        /* Create the FPGA-based scan matcher */
        auto pScanMatcherLikelihood = CreateLikelihoodGreedyEndpoint(
            jsonSettings);
        pScanMatcher = CreateScanMatcherCorrelativeFPGA(
            jsonSettings, std::move(pScanMatcherLikelihood));
    } else {
        /* Create the software scan matcher */
        const std::string scanMatcherType = jsonSettings.get<std::string>(
            "GridMapBuilder.ScanMatcherType");
        auto pScanMatcherLikelihood = CreateLikelihoodGreedyEndpoint(
            jsonSettings);
        pScanMatcher = CreateScanMatcher(
            jsonSettings, scanMatcherType,
            std::move(pScanMatcherLikelihood));
    }

    /* Create the final scan matcher */
    auto pFinalScanMatcherLikelihood = CreateLikelihoodGreedyEndpoint(
        jsonSettings);
    auto pFinalScanMatcher = CreateFinalScanMatcherHillClimbing(
        jsonSettings, std::move(pFinalScanMatcherLikelihood));

    /* Load settings for grid map builder */
    const int numOfParticles = jsonSettings.get<int>(
        "GridMapBuilder.NumOfParticles");

    const double initialPoseX = jsonSettings.get<double>(
        "GridMapBuilder.InitialPose.X");
    const double initialPoseY = jsonSettings.get<double>(
        "GridMapBuilder.InitialPose.Y");
    const double initialPoseTheta = jsonSettings.get<double>(
        "GridMapBuilder.InitialPose.Theta");

    const double mapCellSize = jsonSettings.get<double>(
        "GridMapBuilder.Map.Resolution");
    const bool useLatestMap = jsonSettings.get<bool>(
        "GridMapBuilder.Map.UseLatestMap");
    const std::size_t numOfScansForLatestMap = jsonSettings.get<std::size_t>(
        "GridMapBuilder.Map.NumOfScansForLatestMap");

    const double updateThresholdLinearDist = jsonSettings.get<double>(
        "GridMapBuilder.UpdateThresholdLinearDist");
    const double updateThresholdAngularDist = jsonSettings.get<double>(
        "GridMapBuilder.UpdateThresholdAngularDist");
    const double updateThresholdTime = jsonSettings.get<double>(
        "GridMapBuilder.UpdateThresholdTime");

    const double maxUsableRange = jsonSettings.get<double>(
        "GridMapBuilder.MaxUsableRange");
    const double minUsableRange = jsonSettings.get<double>(
        "GridMapBuilder.MinUsableRange");
    const double resampleThreshold = jsonSettings.get<double>(
        "GridMapBuilder.ResampleThreshold");

    const double probabilityHit = jsonSettings.get<double>(
        "GridMapBuilder.ProbabilityHit");
    const double probabilityMiss = jsonSettings.get<double>(
        "GridMapBuilder.ProbabilityMiss");

    const auto weightNormalizerType = jsonSettings.get<std::string>(
        "GridMapBuilder.WeightNormalizer");
    const auto weightNormalizerHardwareType = jsonSettings.get<std::string>(
        "GridMapBuilder.WeightNormalizerHardware");

    const Mapping::WeightNormalizationType normalizationType =
        useHardwareScanMatcher ?
        DetermineWeightNormalizationType(weightNormalizerHardwareType) :
        DetermineWeightNormalizationType(weightNormalizerType);

    /* Create scan interpolator object if necessary */
    const bool useScanInterpolator = jsonSettings.get<bool>(
        "GridMapBuilder.UseScanInterpolator");
    auto pScanInterpolator = useScanInterpolator ?
        CreateScanInterpolator(jsonSettings) : nullptr;

    /* Create pose covariance estimator */
    auto pCovarianceEstimator = CreateCovarianceEstimator(jsonSettings);

    const double degenerationThreshold = jsonSettings.get<double>(
        "GridMapBuilder.DegenerationThreshold");

    const RobotPose2D<double> initialPose {
        initialPoseX, initialPoseY, initialPoseTheta };

    /* Create the map builder */
    Mapping::MapBuilder mapBuilder {
        maxUsableRange, minUsableRange,
        probabilityHit, probabilityMiss };

    /* Create the grid map builder */
    auto pMapBuilder = std::make_shared<Mapping::GridMapBuilder>(
        std::move(pMotionModel), std::move(pLikelihoodFunc),
        std::move(pScanMatcher), std::move(pFinalScanMatcher),
        std::move(pScanInterpolator), std::move(pCovarianceEstimator),
        std::move(mapBuilder), normalizationType,
        numOfParticles, useLatestMap, numOfScansForLatestMap,
        initialPose, mapCellSize,
        updateThresholdLinearDist, updateThresholdAngularDist,
        updateThresholdTime, resampleThreshold, degenerationThreshold);

    return pMapBuilder;
}

/* Save the metrics */
void SaveMetrics(const std::string& fileName)
{
    /* Save metrics in JSON format */
    auto* const pMetricManager = Metric::MetricManager::Instance();
    /* Convert all metrics to property tree */
    const pt::ptree metricTree = pMetricManager->ToPropertyTree();

    /* Write to JSON file */
    const std::string metricFileName = fileName + ".metric.json";
    pt::write_json(metricFileName, metricTree);
}

/* Load Carmen log data */
void LoadCarmenLog(const fs::path& logFilePath,
                   std::vector<Sensor::SensorDataPtr>& logData)
{
    /* Open Carmen log file */
    std::ifstream logFile { logFilePath };

    /* Exit a program if failed */
    if (!logFile) {
        std::cerr << "Failed to open log file: \'" << logFilePath << "\'\n";
        return;
    }

    /* Load Carmen log file */
    IO::Carmen::CarmenLogReader logReader;
    logReader.Load(logFile, logData);
    logFile.close();
}

/* Struct to store GUI (gnuplot) settings */
struct GuiSettings
{
    bool mEnabled;
    int  mDrawFrameInterval;
};

/* Load the bitstream file to enable the hardware acceleration */
bool LoadBitstream(const pt::ptree& jsonSettings)
{
    /* Load settings for the hardware acceleration */
    const bool enableHardwareAcceleration =
        jsonSettings.get<bool>("Hardware.EnableHardwareAcceleration");
    const std::string bitstreamFileName =
        jsonSettings.get<std::string>("Hardware.BitstreamFileName");

    if (!enableHardwareAcceleration)
        return true;

    /* Load the shared object library for the CMA memory management */
    auto* const pCmaManager = Hardware::CMAMemoryManager::Instance();

    if (!pCmaManager->Load())
        return false;

    /* Load the specified bitstream file */
    Hardware::BitstreamLoader bitstreamLoader;

    if (!bitstreamLoader.Load(bitstreamFileName))
        return false;

    return true;
}

/* Load settings */
bool LoadSettings(const fs::path& settingsFilePath,
                  std::shared_ptr<Mapping::GridMapBuilder>& mapBuilder,
                  GuiSettings& guiSettings)
{
    pt::ptree jsonSettings;
    pt::read_json(settingsFilePath, jsonSettings);

    /* Read settings for GUI */
    guiSettings.mEnabled =
        jsonSettings.get<bool>("GuiEnabled");
    guiSettings.mDrawFrameInterval =
        jsonSettings.get<int>("DrawFrameInterval");

    /* Load the bitstream file to enable the hardware acceleration
     * before setting up the scan matcher */
    if (!LoadBitstream(jsonSettings))
        return false;

    /* Create the grid map builder */
    mapBuilder = CreateGridMapBuilder(jsonSettings);

    return true;
}

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << ' '
                  << "<Carmen log file name>" << ' '
                  << "<Settings file name>" << ' '
                  << "[Output name]" << std::endl;
        return EXIT_FAILURE;
    }

    fs::path settingsFilePath { argv[2] };
    fs::path logFilePath { argv[1] };
    fs::path outputFilePath = (argc == 4) ? argv[3] : logFilePath.stem();

    /* Load Carmen log file */
    std::vector<Sensor::SensorDataPtr> logData;
    LoadCarmenLog(logFilePath, logData);

    if (logData.empty())
        return EXIT_FAILURE;

    /* Load settings from json file */
    std::shared_ptr<Mapping::GridMapBuilder> pMapBuilder;
    GuiSettings guiSettings;

    if (!LoadSettings(settingsFilePath, pMapBuilder, guiSettings))
        return EXIT_FAILURE;

    /* Setup gnuplot helper */
    auto gnuplotHelper = guiSettings.mEnabled ?
        std::make_shared<IO::GnuplotHelper>() : nullptr;

    for (auto& sensorData : logData) {
        auto scanData = std::dynamic_pointer_cast<
            Sensor::ScanData<double>>(std::move(sensorData));

        if (scanData == nullptr)
            continue;

        /* Process scan data */
        const bool mapUpdated = pMapBuilder->ProcessScan(
            scanData, scanData->OdomPose());

        if (!guiSettings.mEnabled || !mapUpdated)
            continue;
        if (pMapBuilder->ProcessCounter() %
            guiSettings.mDrawFrameInterval != 0)
            continue;

        /* Get the trajectory of the best particle */
        const std::size_t bestParticleIdx =
            pMapBuilder->BestParticleIndex();
        const auto bestParticleTrajectory =
            pMapBuilder->ParticleTrajectoryWithTimeStamp(bestParticleIdx);

        /* Draw the trajectory of the best particle */
        gnuplotHelper->DrawParticleTrajectory(bestParticleTrajectory);
    }

    /* Save the map of the best particle */
    const std::size_t bestParticleIdx = pMapBuilder->BestParticleIndex();
    const auto& bestParticle =
        pMapBuilder->ParticleAt(bestParticleIdx);
    const auto bestParticleTrajectory =
        pMapBuilder->ParticleTrajectoryWithTimeStamp(bestParticleIdx);

    IO::MapSaver mapSaver;
    mapSaver.SaveMap(bestParticle.mMap, bestParticleTrajectory,
                     outputFilePath, true);

    /* Save the metrics */
    SaveMetrics(outputFilePath);

    return EXIT_SUCCESS;
}
