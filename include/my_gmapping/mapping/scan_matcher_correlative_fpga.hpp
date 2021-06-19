
/* scan_matcher_correlative_fpga.hpp */

#ifndef MY_GMAPPING_MAPPING_SCAN_MATCHER_CORRELATIVE_FPGA_HPP
#define MY_GMAPPING_MAPPING_SCAN_MATCHER_CORRELATIVE_FPGA_HPP

#include "my_gmapping/hw/ap_ctrl.hpp"
#include "my_gmapping/hw/axi_simple_dma.hpp"
#include "my_gmapping/hw/cma_memory.hpp"
#include "my_gmapping/hw/mmio.hpp"
#include "my_gmapping/mapping/likelihood_function.hpp"
#include "my_gmapping/mapping/scan_matcher.hpp"
#include "my_gmapping/metric/metric.hpp"
#include "my_gmapping/util.hpp"

namespace MyGMapping {
namespace Mapping {

struct ScanMatcherCorrelativeFPGAMetrics
{
    /* `Times` struct holds the processing times for each particle */
    struct Times
    {
        int mInputSetupTime;
        int mSetupIPTime;
        int mScanSendTime;
        int mMapSendTime;
        int mOptimizationTime;
        int mWaitIPTime;
        int mScanMatchingTime;
    };

    /* `Parameters` struct holds the parameter settings for each particle */
    struct Parameters
    {
        /* Metric values related to the algorithmic parameters */
        float mDiffTranslation;
        float mDiffRotation;
        int   mWinSizeX;
        int   mWinSizeY;
        int   mWinSizeTheta;
        float mStepSizeX;
        float mStepSizeY;
        float mStepSizeTheta;
        /* Metric values related to the data transfers */
        bool  mMapTransferred;
        int   mMapChunks;
        bool  mScanTransferred;
        int   mNumOfTransferredScans;
        /* Metric values related to the outputs */
        float mScoreValue;
        float mLikelihoodValue;
    };

    /* Constructor */
    ScanMatcherCorrelativeFPGAMetrics(const std::string& scanMatcherName);
    /* Destructor */
    ~ScanMatcherCorrelativeFPGAMetrics() = default;

    /* Resize the buffer to store the processing times */
    void Resize(const std::size_t numOfParticles);
    /* Set the processing times for each particle */
    void SetTimes(const std::size_t idx, const Times& times);
    /* Set the parameter settings for each particle */
    void SetParameters(const std::size_t idx, const Parameters& params);
    /* Collect the particle-wise metrics and update the overall metrics */
    void Update(const std::size_t bestIdx);

    /* Processing times for all particles */
    std::vector<Times>                mTimes;
    /* Parameter settings for all particles */
    std::vector<Parameters>           mParams;

    /* Total processing time for setting up the input */
    Metric::ValueSequenceBase<int>*   mInputSetupTime;
    /* Total processing time for setting up the scan matcher IP core */
    Metric::ValueSequenceBase<int>*   mSetupIPTime;
    /* Total processing time for sending the scan data */
    Metric::ValueSequenceBase<int>*   mScanSendTime;
    /* Total processing time for sending the grid map data */
    Metric::ValueSequenceBase<int>*   mMapSendTime;
    /* Total processing time for the calculation on the FPGA device */
    Metric::ValueSequenceBase<int>*   mOptimizationTime;
    /* Total processing time for waiting for the scan matcher IP core */
    Metric::ValueSequenceBase<int>*   mWaitIPTime;
    /* Total processing time for the scan matching */
    Metric::ValueSequenceBase<int>*   mScanMatchingTime;

    /* Distance between the initial pose and the final pose */
    Metric::ValueSequenceBase<float>* mDiffTranslation;
    /* Absolute difference between the initial angle and the final angle */
    Metric::ValueSequenceBase<float>* mDiffRotation;
    /* Size of the search window along the X-axis */
    Metric::ValueSequenceBase<int>*   mWinSizeX;
    /* Size of the search window along the Y-axis */
    Metric::ValueSequenceBase<int>*   mWinSizeY;
    /* Size of the search window along the Theta-axis */
    Metric::ValueSequenceBase<int>*   mWinSizeTheta;
    /* Step size along the X-axis */
    Metric::ValueSequenceBase<float>* mStepSizeX;
    /* Step size along the Y-axis */
    Metric::ValueSequenceBase<float>* mStepSizeY;
    /* Step size along the Theta-axis */
    Metric::ValueSequenceBase<float>* mStepSizeTheta;

    /* Number of the grid map transfers */
    Metric::ValueSequenceBase<int>*   mMapTransfers;
    /* Number of the grid map transfer skips */
    Metric::ValueSequenceBase<int>*   mMapTransferSkips;
    /* Number of the transferred grid cell chunks in the grid map */
    Metric::ValueSequenceBase<int>*   mMapChunks;
    /* Number of the scan data transfers */
    Metric::ValueSequenceBase<int>*   mScanTransfers;
    /* Number of the scan data transfer skips */
    Metric::ValueSequenceBase<int>*   mScanTransferSkips;
    /* Number of the transferred scan points */
    Metric::ValueSequenceBase<int>*   mNumOfTransferredScans;

    /* Normalized score value of the best particle */
    Metric::ValueSequenceBase<float>* mScoreValue;
    /* Normalized likelihood value of the best particle */
    Metric::ValueSequenceBase<float>* mLikelihoodValue;
};

/*
 * ScanMatcherCorrelativeFPGAConfig struct holds the necessary information
 * for the real-time correlative scan matcher IP core, namely the base
 * address and address range of the AXI4-Lite slave interface, that are
 * specific for each IP core
 */
struct ScanMatcherCorrelativeFPGAConfig
{
    /* AXI4-Lite slave interface base address */
    std::uint32_t mAxiLiteSBaseAddress;
    /* AXI4-Lite slave interface address range */
    std::uint32_t mAxiLiteSAddressRange;
};

/*
 * ScanMatcherCorrelativeFPGACommonConfig struct holds the necessary
 * information for the real-time correlative scan matcher IP core,
 * including the scan matching parameters, and the address offsets
 * of the control registers, that are shared among all IP cores
 */
struct ScanMatcherCorrelativeFPGACommonConfig
{
    /* Maximum number of the scan points */
    int    mMaxNumOfScans;
    /* Grid map resolution (in meters) */
    double mMapResolution;
    /* Maximum width of the grid map (in the number of grid cells) */
    int    mMaxMapSizeX;
    /* Maximum height of the grid map (in the number of grid cells) */
    int    mMaxMapSizeY;
    /* Resolution of the coarse grid map (in the number of grid cells) */
    int    mLowResolution;
    /* Bit width of the occupancy probability value */
    int    mMapBitWidth;
    /* Width of the grid map chunk (consecutive grid map cells) */
    int    mMapChunkWidth;

    /* Register offsets for the AXI4-Lite slave interface */

    /* AXI4-Lite slave interface control signal */
    std::uint32_t mAxiLiteSApCtrl;
    /* AXI4-Lite slave interface global interrupt enable register */
    std::uint32_t mAxiLiteSGIE;
    /* AXI4-Lite slave interface IP interrupt enable register */
    std::uint32_t mAxiLiteSIER;
    /* AXI4-Lite slave interface IP interrupt status register */
    std::uint32_t mAxiLiteSISR;

    /* Register offsets for the scan matcher */

    /* Register offset for the actual number of the scan points */
    std::uint32_t mAxiLiteSNumOfScans;
    /* Register offset for the maximum scan range considered valid */
    std::uint32_t mAxiLiteSScanRangeMax;
    /* Register offset for the score threshold (for loop detection) */
    std::uint32_t mAxiLiteSScoreThreshold;
    /* Register offset for the sensor pose */
    std::uint32_t mAxiLiteSPoseX;
    std::uint32_t mAxiLiteSPoseY;
    std::uint32_t mAxiLiteSPoseTheta;
    /* Register offset for the actual size of the grid map */
    std::uint32_t mAxiLiteSMapSizeX;
    std::uint32_t mAxiLiteSMapSizeY;
    /* Register offset for the minimum coordinate of the grid map */
    std::uint32_t mAxiLiteSMapMinX;
    std::uint32_t mAxiLiteSMapMinY;
    /* Register offset for the size of the search window */
    std::uint32_t mAxiLiteSWinX;
    std::uint32_t mAxiLiteSWinY;
    std::uint32_t mAxiLiteSWinTheta;
    /* Register offset for the search step */
    std::uint32_t mAxiLiteSStepX;
    std::uint32_t mAxiLiteSStepY;
    std::uint32_t mAxiLiteSStepTheta;
};

/*
 * AxiDmaConfig struct holds the necessary information for the AXI DMA IP core,
 * especially the base address and the address range
 */
struct AxiDmaConfig
{
    /* AXI DMA base address */
    std::uint32_t mBaseAddress;
    /* AXI DMA address range */
    std::uint32_t mAddressRange;
};

class ScanMatcherCorrelativeFPGA final : public ScanMatcher
{
public:
    /* Type definitions */
    using LikelihoodFuncPtr = std::unique_ptr<LikelihoodFunction>;
    using MemoryMappedIOPtr = std::unique_ptr<Hardware::MemoryMappedIO>;
    using AxiSimpleDMAPtr = std::unique_ptr<Hardware::AxiSimpleDMA>;
    using CMAMemory = Hardware::CMAMemory;
    using IPConfig = ScanMatcherCorrelativeFPGAConfig;
    using IPCommonConfig = ScanMatcherCorrelativeFPGACommonConfig;
    using ScanMatcherMetrics = ScanMatcherCorrelativeFPGAMetrics;

    /* Constructor */
    ScanMatcherCorrelativeFPGA(
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
        AxiDmaConfig&& axiDmaConfig1);

    /* Destructor */
    ~ScanMatcherCorrelativeFPGA() = default;

    /* Optimize the particle poses based on the correlative scan matching */
    ScanMatchingResultVector OptimizePose(
        const ScanMatchingQueryVector& queries,
        const Sensor::ScanDataPtr<double>& scanData) override;

private:
    /* Optimize the particle poses using the scan matcher IP core */
    ScanMatchingResultVector OptimizePoseCore(
        const int coreId,
        const std::size_t idxBegin,
        const std::size_t idxEnd,
        const ScanMatchingQueryVector& queries,
        const Sensor::ScanDataPtr<double>& scanData);

    /* Compute the search step */
    void ComputeSearchStep(const double mapResolution,
                           const Sensor::ScanDataPtr<double>& scanData,
                           double& stepX,
                           double& stepY,
                           double& stepTheta) const;

    /* Initialize the scan matcher IP core */
    void Initialize(const int coreId);
    /* Initialize the CMA memory for the input data */
    void InitializeCMAMemoryInput(const int coreId);
    /* Initialize the CMA memory for the output data */
    void InitializeCMAMemoryOutput(const int coreId);

    /* Compute the bounding box of the grid map */
    BoundingBox<int> ComputeBoundingBox(
        const GridMap& gridMap,
        const RobotPose2D<double>& sensorPose) const;

    /* Set the scan matching parameters through AXI4-Lite slave interface */
    void SetParameterRegisters(
        const int coreId,
        const int numOfScans,
        const double scanRangeMax,
        const int scoreThreshold,
        const RobotPose2D<double>& minSensorPose,
        const Point2D<int>& gridMapSize,
        const Point2D<double>& gridMapMinPos,
        const int winX, const int winY, const int winTheta,
        const double stepX, const double stepY, const double stepTheta);

    /* Send the scan data through AXI DMA */
    void SendScanData(const int coreId,
                      const bool scanDataTransferred,
                      const Sensor::ScanDataPtr<double>& scanData,
                      ScanMatcherMetrics::Parameters& params);
    /* Send the grid map through AXI DMA */
    void SendGridMap(const int coreId,
                     const GridMap& gridMap,
                     const BoundingBox<int>& desiredBox,
                     ScanMatcherMetrics::Parameters& params);
    /* Receive the result through AXI DMA */
    void ReceiveResult(const int coreId,
                       int& scoreMax, int& bestX, int& bestY, int& bestTheta);

    /* Start the scan matcher IP core */
    void StartIPCore(const int coreId);
    /* Wait until the scan matcher IP core is in idle state */
    void WaitIPCore(const int coreId);

    /* Read the control register of the scan matcher IP core */
    inline std::uint32_t ReadCtrlReg(const int coreId)
    { return this->mControlRegisters[coreId]->Read(
        this->mCommonConfig.mAxiLiteSApCtrl); }
    /* Write to the control register of the scan matcher IP core */
    inline void WriteCtrlReg(const int coreId, std::uint32_t value)
    { this->mControlRegisters[coreId]->Write(
        this->mCommonConfig.mAxiLiteSApCtrl, value); }

private:
    /* Likelihood function to compute the observation likelihood */
    const LikelihoodFuncPtr mLikelihoodFunc;
    /* Linear (horizontal) size of the searching window */
    const double            mRangeX;
    /* Linear (vertical) size of the searching window */
    const double            mRangeY;
    /* Angular range of the searching window */
    const double            mRangeTheta;
    /* Maximum laser scan range considered for scan matching */
    const double            mScanRangeMax;

    /* Total number of the real-time correlative-based scan matcher
     * IP cores implemented on the Zynq device */
    static constexpr const int        NumOfIPCores = 2;

    /* Configuration of the scan matcher IP core */
    const IPConfig       mConfig[NumOfIPCores];
    const IPCommonConfig mCommonConfig;
    /* Configuration of the AXI DMA IP core */
    const AxiDmaConfig   mAxiDmaConfig[NumOfIPCores];
    /* Interface to the scan matcher IP core control registers */
    MemoryMappedIOPtr    mControlRegisters[NumOfIPCores];
    /* Interface to the AXI DMA IP core */
    AxiSimpleDMAPtr      mAxiDma[NumOfIPCores];
    /* CMA memory to store the scan matching input */
    CMAMemory            mInputData[NumOfIPCores];
    /* CMA memory to store the scan matching result */
    CMAMemory            mOutputData[NumOfIPCores];
    /* Metrics information */
    ScanMatcherMetrics   mMetrics;
};

} /* namespace Mapping */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_MAPPING_SCAN_MATCHER_CORRELATIVE_FPGA_HPP */
