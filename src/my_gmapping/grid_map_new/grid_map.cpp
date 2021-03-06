
/* grid_map.cpp */

#include "my_gmapping/grid_map_new/grid_map.hpp"
#include "my_gmapping/grid_map_new/grid_binary_bayes.hpp"
#include "my_gmapping/grid_map_new/grid_constant.hpp"
#include "my_gmapping/grid_map_new/grid_counted.hpp"

#if (defined __ARM_NEON) || (defined __ARM_NEON__)
#include <arm_neon.h>
#endif

namespace MyGMapping {
namespace GridMapNew {

/* Template class declarations */
template class GridMap<GridBinaryBayes>;
template class GridMap<GridConstant>;
template class GridMap<GridCounted>;

/* Round up to the nearest power of 2
 * Borrowed from https://stackoverflow.com/questions/466204 */
int ToNearestPowerOf2(int x)
{
    Assert(x > 0);

    --x;
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);
    ++x;

    return x;
}

/* Constructor to create an empty grid map */
template <typename T>
GridMap<T>::GridMap() :
    BaseType(),
    mLog2BlockSize(0),
    mBlockSize(0),
    mBlockRows(0),
    mBlockCols(0),
    mBlocks(nullptr)
{
}

/* Constructor with the number of rows and columns */
template <typename T>
GridMap<T>::GridMap(const double resolution, const int blockSize,
                    const int desiredRows, const int desiredCols) :
    BaseType(),
    mLog2BlockSize(0),
    mBlockSize(0),
    mBlockRows(0),
    mBlockCols(0),
    mBlocks(nullptr)
{
    Assert(resolution > 0.0);
    Assert(blockSize > 0);
    Assert(desiredRows > 0);
    Assert(desiredCols > 0);

    this->Initialize(resolution, blockSize, desiredRows, desiredCols);
}

/* Constructor with the width and height */
template <typename T>
GridMap<T>::GridMap(const double resolution, const int blockSize,
                    const double width, const double height) :
    BaseType(),
    mLog2BlockSize(0),
    mBlockSize(0),
    mBlockRows(0),
    mBlockCols(0),
    mBlocks(nullptr)
{
    Assert(resolution > 0.0);
    Assert(blockSize > 0);
    Assert(width > 0.0);
    Assert(height > 0.0);

    /* Compute the number of desired rows and columns */
    const int desiredRows = static_cast<int>(std::ceil(height / resolution));
    const int desiredCols = static_cast<int>(std::ceil(width / resolution));

    this->Initialize(resolution, blockSize, desiredRows, desiredCols);
}

/* Constructor with the number of block rows, block columns,
 * and the offset from the current map-local coordinate frame
 * to the original map-local coordinate frame */
template <typename T>
GridMap<T>::GridMap(const double resolution, const int blockSize,
                    const int blockRows, const int blockCols,
                    const Point2D<double>& posOffset) :
    BaseType(),
    mLog2BlockSize(0),
    mBlockSize(0),
    mBlockRows(0),
    mBlockCols(0),
    mBlocks(nullptr)
{
    Assert(blockSize > 0);
    Assert(blockRows > 0);
    Assert(blockCols > 0);

    /* Compute the base-2 logarithm of the block size */
    const int powerOf2 = ToNearestPowerOf2(blockSize);
    const int log2BlockSize = __builtin_ctz(powerOf2);

    /* Set the size of the blocks */
    this->mLog2BlockSize = log2BlockSize;
    this->mBlockSize = blockSize;
    this->mBlockRows = blockRows;
    this->mBlockCols = blockCols;

    /* Allocate the number of blocks */
    this->Allocate();

    /* Initialize the geometric information of the grid map */
    const int rows = this->mBlockRows << this->mLog2BlockSize;
    const int cols = this->mBlockCols << this->mLog2BlockSize;
    this->mGeometry.Initialize(resolution, rows, cols);

    /* Set the positional offset */
    this->mGeometry.mPosOffset = posOffset;
}

/* Copy constructor */
template <typename T>
GridMap<T>::GridMap(const GridMap& other) :
    BaseType(other),
    mLog2BlockSize(other.mLog2BlockSize),
    mBlockSize(other.mBlockSize),
    mBlockRows(other.mBlockRows),
    mBlockCols(other.mBlockCols),
    mBlocks(nullptr)
{
    if (other.mBlocks == nullptr)
        return;

    /* Allocate the storage for blocks */
    this->Allocate();

    /* Copy the blocks */
    const int numOfBlocks = this->mBlockRows * this->mBlockCols;
    std::copy_n(other.mBlocks.get(), numOfBlocks, this->mBlocks.get());
}

/* Copy assignment operator */
template <typename T>
GridMap<T>& GridMap<T>::operator=(const GridMap& other)
{
    if (this == &other)
        return *this;

    /* Release the storage for blocks if not valid */
    if (other.mBlocks == nullptr) {
        this->Reset();
        return *this;
    }

    /* Reallocate the storage for blocks if not valid */
    if (this->mLog2BlockSize != other.mLog2BlockSize ||
        this->mBlockRows != other.mBlockRows ||
        this->mBlockCols != other.mBlockCols) {
        this->mLog2BlockSize = other.mLog2BlockSize;
        this->mBlockSize = other.mBlockSize;
        this->mBlockRows = other.mBlockRows;
        this->mBlockCols = other.mBlockCols;
        this->Allocate();
    }

    /* Copy the blocks */
    const int numOfBlocks = this->mBlockRows * this->mBlockCols;
    std::copy_n(other.mBlocks.get(), numOfBlocks, this->mBlocks.get());

    /* Copy the geometric information */
    this->mGeometry = other.mGeometry;

    return *this;
}

/* Move constructor */
template <typename T>
GridMap<T>::GridMap(GridMap&& other) noexcept :
    BaseType(std::move(other)),
    mLog2BlockSize(other.mLog2BlockSize),
    mBlockSize(other.mBlockSize),
    mBlockRows(other.mBlockRows),
    mBlockCols(other.mBlockCols),
    mBlocks(std::move(other.mBlocks))
{
}

/* Move assignment operator */
template <typename T>
GridMap<T>& GridMap<T>::operator=(GridMap&& other) noexcept
{
    if (this == &other)
        return *this;

    BaseType::operator=(std::move(other));
    this->mLog2BlockSize = other.mLog2BlockSize;
    this->mBlockSize = other.mBlockSize;
    this->mBlockRows = other.mBlockRows;
    this->mBlockCols = other.mBlockCols;
    this->mBlocks = std::move(other.mBlocks);

    return *this;
}

/* Initialize with the number of rows and columns */
template <typename T>
void GridMap<T>::Initialize(const double resolution, const int blockSize,
                            const int desiredRows, const int desiredCols)
{
    /* Compute the base-2 logarithm of the block size */
    const int powerOf2 = ToNearestPowerOf2(blockSize);
    const int log2BlockSize = __builtin_ctz(powerOf2);

    /* Set the size of the blocks */
    this->mLog2BlockSize = log2BlockSize;
    this->mBlockSize = 1 << log2BlockSize;

    /* Compute the number of blocks */
    this->mBlockRows = (desiredRows + this->mBlockSize - 1)
                       >> this->mLog2BlockSize;
    this->mBlockCols = (desiredCols + this->mBlockSize - 1)
                       >> this->mLog2BlockSize;

    /* Allocate the storage for the blocks */
    this->Allocate();

    /* Initialize the geometric information of the grid map */
    const int rows = this->mBlockRows << this->mLog2BlockSize;
    const int cols = this->mBlockCols << this->mLog2BlockSize;
    this->mGeometry.Initialize(resolution, rows, cols);
}

/* Initialize with the size and the position offset */
template <typename T>
void GridMap<T>::Initialize(const double resolution, const int blockSize,
                            const int desiredRows, const int desiredCols,
                            const double offsetX, const double offsetY)
{
    /* Compute the base-2 logarithm of the block size */
    const int powerOf2 = ToNearestPowerOf2(blockSize);
    const int log2BlockSize = __builtin_ctz(powerOf2);

    /* Set the size of the blocks */
    this->mLog2BlockSize = log2BlockSize;
    this->mBlockSize = 1 << log2BlockSize;

    /* Compute the number of blocks */
    this->mBlockRows = (desiredRows + this->mBlockSize - 1)
                       >> this->mLog2BlockSize;
    this->mBlockCols = (desiredCols + this->mBlockSize - 1)
                       >> this->mLog2BlockSize;

    /* Allocate the storage for the blocks */
    this->Allocate();

    /* Initialize the geometric information of the grid map */
    const int rows = this->mBlockRows << this->mLog2BlockSize;
    const int cols = this->mBlockCols << this->mLog2BlockSize;
    this->mGeometry.Initialize(resolution, rows, cols, offsetX, offsetY);
}

/* Reset to the initial state */
template <typename T>
void GridMap<T>::Reset()
{
    this->mGeometry.Reset();
    this->mLog2BlockSize = 0;
    this->mBlockSize = 0;
    this->mBlockRows = 0;
    this->mBlockCols = 0;
    this->Release();
}

/* Allocate the storage for the blocks */
template <typename T>
void GridMap<T>::Allocate()
{
    const int numOfBlocks = this->mBlockRows * this->mBlockCols;
    this->mBlocks.reset(new T[numOfBlocks]);
    Assert(this->mBlocks != nullptr);
}

/* Release the storage for the blocks */
template <typename T>
void GridMap<T>::Release()
{
    this->mBlocks.reset(nullptr);
}

/* Reset the grid values to unknown */
template <typename T>
void GridMap<T>::ResetValues()
{
    const int numOfBlocks = this->mBlockRows * this->mBlockCols;
    T* block = this->mBlocks.get();

    for (int i = 0; i < numOfBlocks; ++i, ++block)
        if (block->IsAllocated())
            block->ResetValues();
}

/* Copy the internal values to the buffer */
template <typename T>
template <typename U, typename F>
void GridMap<T>::CopyValuesInternal(
    U* buffer, const BoundingBox<int>& boundingBox, F copyFunc) const
{
    Assert(boundingBox.mMin.mX >= 0);
    Assert(boundingBox.mMin.mY >= 0);
    Assert(boundingBox.mMax.mX <= this->mGeometry.mCols);
    Assert(boundingBox.mMax.mY <= this->mGeometry.mRows);
    Assert(boundingBox.mMax.mX > boundingBox.mMin.mX);
    Assert(boundingBox.mMax.mY > boundingBox.mMin.mY);

    Assert(this->mBlockSize % sizeof(U) == 0);
    Assert(boundingBox.mMin.mX % sizeof(U) == 0);
    Assert(boundingBox.mMax.mX % sizeof(U) == 0);

    /* Buffer should have the same size as the bounding box */

    /* Convert the grid index to the block index */
    const Point2D<int> blockIdxMin =
        this->IndexToBlock(boundingBox.mMin.mY, boundingBox.mMin.mX);
    const Point2D<int> blockIdxMax =
        this->IndexToBlock(boundingBox.mMax.mY + this->mBlockSize - 1,
                           boundingBox.mMax.mX + this->mBlockSize - 1);

    const int cols = boundingBox.mMax.mX - boundingBox.mMin.mX;
    const T* block = this->Block(blockIdxMin.mY, blockIdxMin.mX);

    for (int row = blockIdxMin.mY; row < blockIdxMax.mY; ++row) {
        const int blockRowMin = (row == blockIdxMin.mY) ?
            boundingBox.mMin.mY - (row << this->mLog2BlockSize) : 0;
        const int blockRowMax = (row == blockIdxMax.mY - 1) ?
            boundingBox.mMax.mY - (row << this->mLog2BlockSize) :
            this->mBlockSize;

        auto* currentBuffer = buffer;
        const auto* currentBlock = block;

        for (int col = blockIdxMin.mX; col < blockIdxMax.mX; ++col) {
            const int blockColMin = (col == blockIdxMin.mX) ?
                boundingBox.mMin.mX - (col << this->mLog2BlockSize) : 0;
            const int blockColMax = (col == blockIdxMax.mX - 1) ?
                boundingBox.mMax.mX - (col << this->mLog2BlockSize) :
                this->mBlockSize;
            const BoundingBox<int> blockBoundingBox {
                blockColMin, blockRowMin, blockColMax, blockRowMax };

            /* Copy the grid values to the buffer */
            copyFunc(currentBlock, currentBuffer, cols, blockBoundingBox);

            /* Advance the pointer to the block */
            currentBlock++;
            /* Advance the pointer to the buffer */
            currentBuffer += ((blockColMax - blockColMin) / sizeof(U));
        }

        /* Advance the pointer to the block */
        block += this->mBlockCols;
        /* Advance the pointer to the buffer */
        buffer += (cols * (blockRowMax - blockRowMin) / sizeof(U));
    }
}

/* Check if the grid cell is allocated on the heap */
template <typename T>
bool GridMap<T>::IsAllocated(const int row, const int col) const
{
    const Point2D<int> blockIdx = this->IndexToBlock(row, col);
    return this->IsBlockAllocated(blockIdx.mY, blockIdx.mX);
}

/* Get the internal value of the grid cell */
template <typename T>
typename T::ValueType GridMap<T>::Value(
    const int row, const int col) const
{
    const Point2D<int> blockIdx = this->IndexToBlock(row, col);
    Assert(this->IsBlockAllocated(blockIdx.mY, blockIdx.mX));

    const Point2D<int> blockOffset = this->IndexToBlockOffset(row, col);
    const T* block = this->Block(blockIdx.mY, blockIdx.mX);
    return block->ValueUnchecked(blockOffset.mY, blockOffset.mX);
}

/* Get the internal value of the grid cell (without checks) */
template <typename T>
typename T::ValueType GridMap<T>::ValueUnchecked(
    const int row, const int col) const
{
    const Point2D<int> blockIdx = this->IndexToBlock(row, col);
    const Point2D<int> blockOffset = this->IndexToBlockOffset(row, col);
    const T* block = this->Block(blockIdx.mY, blockIdx.mX);
    return block->ValueUnchecked(blockOffset.mY, blockOffset.mX);
}

/* Get the internal value of the grid cell or return the default value */
template <typename T>
typename T::ValueType GridMap<T>::ValueOr(
    const int row, const int col, const typename T::ValueType value) const
{
    const Point2D<int> blockIdx = this->IndexToBlock(row, col);

    if (!this->IsBlockAllocated(blockIdx.mY, blockIdx.mX))
        return value;

    const Point2D<int> blockOffset = this->IndexToBlockOffset(row, col);
    const T* block = this->Block(blockIdx.mY, blockIdx.mX);
    return block->ValueUnchecked(blockOffset.mY, blockOffset.mX);
}

/* Get the probability value of the grid cell */
template <typename T>
typename T::ProbabilityType GridMap<T>::Probability(
    const int row, const int col) const
{
    const Point2D<int> blockIdx = this->IndexToBlock(row, col);
    Assert(this->IsBlockAllocated(blockIdx.mY, blockIdx.mX));

    const Point2D<int> blockOffset = this->IndexToBlockOffset(row, col);
    const T* block = this->Block(blockIdx.mY, blockIdx.mX);
    return block->ProbabilityUnchecked(blockOffset.mY, blockOffset.mX);
}

/* Get the probability value of the grid cell (without checks) */
template <typename T>
typename T::ProbabilityType GridMap<T>::ProbabilityUnchecked(
    const int row, const int col) const
{
    const Point2D<int> blockIdx = this->IndexToBlock(row, col);
    const Point2D<int> blockOffset = this->IndexToBlockOffset(row, col);
    const T* block = this->Block(blockIdx.mY, blockIdx.mX);
    return block->ProbabilityUnchecked(blockOffset.mY, blockOffset.mX);
}

/* Get the probability value of the grid cell or return the default value */
template <typename T>
typename T::ProbabilityType GridMap<T>::ProbabilityOr(
    const int row, const int col, const typename T::ProbabilityType prob) const
{
    const Point2D<int> blockIdx = this->IndexToBlock(row, col);

    if (!this->IsBlockAllocated(blockIdx.mY, blockIdx.mX))
        return prob;

    const Point2D<int> blockOffset = this->IndexToBlockOffset(row, col);
    const T* block = this->Block(blockIdx.mY, blockIdx.mX);
    return block->ProbabilityUnchecked(blockOffset.mY, blockOffset.mX);
}

/* Copy the internal values to the buffer */
template <typename T>
void GridMap<T>::CopyValues(typename T::ValueType* buffer) const
{
    const BoundingBox<int> boundingBox {
        0, 0, this->mGeometry.mCols, this->mGeometry.mRows };
    this->CopyValues(buffer, boundingBox);
}

/* Copy the internal values to the buffer */
template <typename T>
void GridMap<T>::CopyValues(typename T::ValueType* buffer,
                            const BoundingBox<int>& boundingBox) const
{
    const auto copyFunc = [](const T* currentBlock,
                             typename T::ValueType* currentBuffer,
                             int bufferCols, const BoundingBox<int>& box) {
        currentBlock->CopyValues(currentBuffer, bufferCols, box); };
    this->CopyValuesInternal(buffer, boundingBox, copyFunc);
}

/* Copy the internal values as std::uint8_t to the buffer */
template <typename T>
void GridMap<T>::CopyValuesU8(std::uint8_t* buffer) const
{
    const BoundingBox<int> boundingBox {
        0, 0, this->mGeometry.mCols, this->mGeometry.mRows };
    this->CopyValuesU8(buffer, boundingBox);
}

/* Copy the internal values in the specified region as std::uint8_t */
template <typename T>
void GridMap<T>::CopyValuesU8(std::uint8_t* buffer,
                              const BoundingBox<int>& boundingBox) const
{
    const auto copyFunc = [](const T* currentBlock,
                             std::uint8_t* currentBuffer,
                             int bufferCols, const BoundingBox<int>& box) {
        currentBlock->CopyValuesU8(currentBuffer, bufferCols, box); };
    this->CopyValuesInternal(buffer, boundingBox, copyFunc);
}

/* Copy the internal values as std::uint8_t to the buffer
 * and ensure the 4-byte aligned accesses */
template <typename T>
void GridMap<T>::CopyValuesU8x4(std::uint32_t* buffer) const
{
    const BoundingBox<int> boundingBox {
        0, 0, this->mGeometry.mCols, this->mGeometry.mRows };
    this->CopyValuesU8x4(buffer, boundingBox);
}

/* Copy the internal values as std::uint8_t to the buffer
 * and ensure the 4-byte aligned accesses */
template <typename T>
void GridMap<T>::CopyValuesU8x4(std::uint32_t* buffer,
                                const BoundingBox<int>& boundingBox) const
{
    const auto copyFunc = [](const T* currentBlock,
                             std::uint32_t* currentBuffer,
                             int bufferCols, const BoundingBox<int>& box) {
        currentBlock->CopyValuesU8x4(currentBuffer, bufferCols, box); };
    this->CopyValuesInternal(buffer, boundingBox, copyFunc);
}

/* Set the internal value of the grid cell */
template <typename T>
void GridMap<T>::SetValue(
    const int row, const int col, const typename T::ValueType value)
{
    const Point2D<int> blockIdx = this->IndexToBlock(row, col);
    Assert(this->IsBlockInside(blockIdx.mY, blockIdx.mX));

    const Point2D<int> blockOffset = this->IndexToBlockOffset(row, col);
    T* block = this->Block(blockIdx.mY, blockIdx.mX);

    /* Allocate and initialize the block if necessary */
    if (!block->IsAllocated())
        block->Initialize(this->mLog2BlockSize);

    block->SetValueUnchecked(blockOffset.mY, blockOffset.mX, value);
}

/* Set the internal value of the grid cell (without checks) */
template <typename T>
void GridMap<T>::SetValueUnchecked(
    const int row, const int col, const typename T::ValueType value)
{
    const Point2D<int> blockIdx = this->IndexToBlock(row, col);
    const Point2D<int> blockOffset = this->IndexToBlockOffset(row, col);
    T* block = this->Block(blockIdx.mY, blockIdx.mX);

    /* Allocate and initialize the block if necessary */
    if (!block->IsAllocated())
        block->Initialize(this->mLog2BlockSize);

    block->SetValueUnchecked(blockOffset.mY, blockOffset.mX, value);
}

/* Set the probability value of the grid cell */
template <typename T>
void GridMap<T>::SetProbability(
    const int row, const int col, const typename T::ProbabilityType prob)
{
    const Point2D<int> blockIdx = this->IndexToBlock(row, col);
    Assert(this->IsBlockInside(blockIdx.mY, blockIdx.mX));

    const Point2D<int> blockOffset = this->IndexToBlockOffset(row, col);
    T* block = this->Block(blockIdx.mY, blockIdx.mX);

    /* Allocate and initialize the block if necessary */
    if (!block->IsAllocated())
        block->Initialize(this->mLog2BlockSize);

    block->SetProbabilityUnchecked(blockOffset.mY, blockOffset.mX, prob);
}

/* Set the probability value of the grid cell (without checks) */
template <typename T>
void GridMap<T>::SetProbabilityUnchecked(
    const int row, const int col, const typename T::ProbabilityType prob)
{
    const Point2D<int> blockIdx = this->IndexToBlock(row, col);
    const Point2D<int> blockOffset = this->IndexToBlockOffset(row, col);
    T* block = this->Block(blockIdx.mY, blockIdx.mX);

    /* Allocate and initialize the block if necessary */
    if (!block->IsAllocated())
        block->Initialize(this->mLog2BlockSize);

    block->SetProbabilityUnchecked(blockOffset.mY, blockOffset.mX, prob);
}

/* Fill all grid values with the given internal value */
template <typename T>
void GridMap<T>::FillValue(const typename T::ValueType value)
{
    const int numOfBlocks = this->mBlockRows * this->mBlockCols;
    T* block = this->mBlocks.get();

    for (int i = 0; i < numOfBlocks; ++i, ++block)
        if (block->IsAllocated())
            block->FillValue(value);
}

/* Fill all grid values with the given probability value */
template <typename T>
void GridMap<T>::FillProbability(const typename T::ProbabilityType prob)
{
    const int numOfBlocks = this->mBlockRows * this->mBlockCols;
    T* block = this->mBlocks.get();

    for (int i = 0; i < numOfBlocks; ++i, ++block)
        if (block->IsAllocated())
            block->FillProbability(prob);
}

/* Update the grid value given an observation */
template <typename T>
void GridMap<T>::Update(const int row, const int col,
                        const typename T::ObservationType observation)
{
    const Point2D<int> blockIdx = this->IndexToBlock(row, col);
    Assert(this->IsBlockInside(blockIdx.mY, blockIdx.mX));

    const Point2D<int> blockOffset = this->IndexToBlockOffset(row, col);
    T* block = this->Block(blockIdx.mY, blockIdx.mX);

    /* Allocate and initialize the block if necessary */
    if (!block->IsAllocated())
        block->Initialize(this->mLog2BlockSize);

    block->UpdateUnchecked(blockOffset.mY, blockOffset.mX, observation);
}

/* Update the grid value given an observation (without input checks) */
template <typename T>
void GridMap<T>::UpdateUnchecked(const int row, const int col,
                                 const typename T::ObservationType observation)
{
    const Point2D<int> blockIdx = this->IndexToBlock(row, col);
    const Point2D<int> blockOffset = this->IndexToBlockOffset(row, col);
    T* block = this->Block(blockIdx.mY, blockIdx.mX);

    /* Allocate and initialize the block if necessary */
    if (!block->IsAllocated())
        block->Initialize(this->mLog2BlockSize);

    block->UpdateUnchecked(blockOffset.mY, blockOffset.mX, observation);
}

/* Update the grid value given an odds */
template <typename T>
void GridMap<T>::UpdateOdds(const int row, const int col,
                            const typename T::ObservationType odds)
{
    const Point2D<int> blockIdx = this->IndexToBlock(row, col);
    Assert(this->IsBlockInside(blockIdx.mY, blockIdx.mX));

    const Point2D<int> blockOffset = this->IndexToBlockOffset(row, col);
    T* block = this->Block(blockIdx.mY, blockIdx.mX);

    /* Allocate and initialize the block if necessary */
    if (!block->IsAllocated())
        block->Initialize(this->mLog2BlockSize);

    block->UpdateOddsUnchecked(blockOffset.mY, blockOffset.mX, odds);
}

/* Update the grid value given an odds (without input checks) */
template <typename T>
void GridMap<T>::UpdateOddsUnchecked(const int row, const int col,
                                     const typename T::ObservationType odds)
{
    const Point2D<int> blockIdx = this->IndexToBlock(row, col);
    const Point2D<int> blockOffset = this->IndexToBlockOffset(row, col);
    T* block = this->Block(blockIdx.mY, blockIdx.mX);

    /* Allocate and initialize the block if necessary */
    if (!block->IsAllocated())
        block->Initialize(this->mLog2BlockSize);

    block->UpdateOddsUnchecked(blockOffset.mY, blockOffset.mX, odds);
}

/* Update multiple grid values given an odds */
template <typename T>
void GridMap<T>::UpdateOdds(const std::vector<Point2D<int>>& indices,
                            const typename T::ObservationType odds)
{
    /* Update each grid cell value using the same odds */
    for (const auto& idx : indices) {
        const auto blockIdx = this->IndexToBlock(idx.mY, idx.mX);
        Assert(this->IsBlockInside(blockIdx.mY, blockIdx.mX));

        const auto blockOffset = this->IndexToBlockOffset(idx.mY, idx.mX);
        T* block = this->Block(blockIdx.mY, blockIdx.mX);

        /* Allocate and initialize the block if necessary */
        if (!block->IsAllocated())
            block->Initialize(this->mLog2BlockSize);

        block->UpdateOddsUnchecked(blockOffset.mY, blockOffset.mX, odds);
    }
}

#if (defined __ARM_NEON) || (defined __ARM_NEON__)

/* Update multiple grid values given an odds (without input checks) */
template <typename T>
void GridMap<T>::UpdateOddsUnchecked(const std::vector<Point2D<int>>& indices,
                                     const typename T::ObservationType odds)
{
    const std::size_t sizeRounded = (indices.size() >> 2) << 2;
    const int mask = (1 << this->mLog2BlockSize) - 1;
    const int32x4_t mask4 = vmovq_n_s32(mask);
    const int32x4_t shift4 = vmovq_n_s32(-this->mLog2BlockSize);

    /* Update each grid cell value using the same odds */
    for (std::size_t i = 0; i < sizeRounded; i += 4) {
        /* Assume that these indices are valid, i.e.,
         * 0 <= col < this->mGeometry.mCols and
         * 0 <= row < this->mGeometry.mRows */
        const int32x4_t col4 { indices[i].mX, indices[i + 1].mX,
                               indices[i + 2].mX, indices[i + 3].mX };
        const int32x4_t row4 { indices[i].mY, indices[i + 1].mY,
                               indices[i + 2].mY, indices[i + 3].mY };

        /* Compute the block indices and offsets */
        /* These operations do not care about the negative indices */
        const int32x4_t blockRow4 = vshlq_s32(row4, shift4);
        const int32x4_t blockCol4 = vshlq_s32(col4, shift4);
        const int32x4_t blockRowOffset4 = vandq_s32(row4, mask4);
        const int32x4_t blockColOffset4 = vandq_s32(col4, mask4);

        T* blocks[4] = { this->Block(vgetq_lane_s32(blockRow4, 0),
                                     vgetq_lane_s32(blockCol4, 0)),
                         this->Block(vgetq_lane_s32(blockRow4, 1),
                                     vgetq_lane_s32(blockCol4, 1)),
                         this->Block(vgetq_lane_s32(blockRow4, 2),
                                     vgetq_lane_s32(blockCol4, 2)),
                         this->Block(vgetq_lane_s32(blockRow4, 3),
                                     vgetq_lane_s32(blockCol4, 3)) };

        /* Allocate and initialize the block if necessary */
        if (!blocks[0]->IsAllocated())
            blocks[0]->Initialize(this->mLog2BlockSize);
        if (!blocks[1]->IsAllocated())
            blocks[1]->Initialize(this->mLog2BlockSize);
        if (!blocks[2]->IsAllocated())
            blocks[2]->Initialize(this->mLog2BlockSize);
        if (!blocks[3]->IsAllocated())
            blocks[3]->Initialize(this->mLog2BlockSize);

        blocks[0]->UpdateOddsUnchecked(vgetq_lane_s32(blockRowOffset4, 0),
                                       vgetq_lane_s32(blockColOffset4, 0),
                                       odds);
        blocks[1]->UpdateOddsUnchecked(vgetq_lane_s32(blockRowOffset4, 1),
                                       vgetq_lane_s32(blockColOffset4, 1),
                                       odds);
        blocks[2]->UpdateOddsUnchecked(vgetq_lane_s32(blockRowOffset4, 2),
                                       vgetq_lane_s32(blockColOffset4, 2),
                                       odds);
        blocks[3]->UpdateOddsUnchecked(vgetq_lane_s32(blockRowOffset4, 3),
                                       vgetq_lane_s32(blockColOffset4, 3),
                                       odds);
    }

    for (std::size_t i = sizeRounded; i < indices.size(); ++i) {
        const int row = indices[i].mY;
        const int col = indices[i].mX;
        const int blockRow = row >> this->mLog2BlockSize;
        const int blockCol = col >> this->mLog2BlockSize;
        const int blockRowOffset = row & mask;
        const int blockColOffset = col & mask;

        T* block = this->Block(blockRow, blockCol);

        /* Allocate and initialize the block if necessary */
        if (!block->IsAllocated())
            block->Initialize(this->mLog2BlockSize);

        block->UpdateOddsUnchecked(blockRowOffset, blockColOffset, odds);
    }
}

#else

/* Update multiple grid values given an odds (without input checks) */
template <typename T>
void GridMap<T>::UpdateOddsUnchecked(const std::vector<Point2D<int>>& indices,
                                     const typename T::ObservationType odds)
{
    /* Update each grid cell value using the same odds */
    for (const auto& idx : indices) {
        const auto blockIdx = this->IndexToBlock(idx.mY, idx.mX);
        const auto blockOffset = this->IndexToBlockOffset(idx.mY, idx.mX);
        T* block = this->Block(blockIdx.mY, blockIdx.mX);

        /* Allocate and initialize the block if necessary */
        if (!block->IsAllocated())
            block->Initialize(this->mLog2BlockSize);

        block->UpdateOddsUnchecked(blockOffset.mY, blockOffset.mX, odds);
    }
}

#endif

/* Check if the block index is valid */
template <typename T>
bool GridMap<T>::IsBlockInside(const int row, const int col) const
{
    return (row >= 0 && row < this->mBlockRows) &&
           (col >= 0 && col < this->mBlockCols);
}

/* Check if the block is allocated */
template <typename T>
bool GridMap<T>::IsBlockAllocated(const int row, const int col) const
{
    return this->IsBlockInside(row, col) &&
           this->Block(row, col)->IsAllocated();
}

/* Convert the grid index to the block index */
template <typename T>
Point2D<int> GridMap<T>::IndexToBlock(const int row, const int col) const
{
    /* The below code assumes that the right-shift on the signed integer is
     * the arithmetic right shift and not the logical (unsigned) shift */
    const int blockRow = (row >= 0) ? (row >> this->mLog2BlockSize) :
                         ((row >> this->mLog2BlockSize) - 1);
    const int blockCol = (col >= 0) ? (col >> this->mLog2BlockSize) :
                         ((col >> this->mLog2BlockSize) - 1);
    return Point2D<int> { blockCol, blockRow };
}

/* Convert the grid index to the block index offset */
template <typename T>
Point2D<int> GridMap<T>::IndexToBlockOffset(const int row, const int col) const
{
    /* The below code assumes that the negative numbers are expressed in
     * the two's complement (the offset for -1 is 31 and for -2 is 30) and
     * the block size is not too large */
    const int mask = (1 << this->mLog2BlockSize) - 1;
    const int offsetRow = row & mask;
    const int offsetCol = col & mask;
    return Point2D<int> { offsetCol, offsetRow };
}

/* Convert the block index to the grid index range */
template <typename T>
BoundingBox<int> GridMap<T>::BlockToIndexRange(
    const int row, const int col) const
{
    const int rowMin = row << this->mLog2BlockSize;
    const int colMin = col << this->mLog2BlockSize;
    const int rowMax = (row + 1) << this->mLog2BlockSize;
    const int colMax = (col + 1) << this->mLog2BlockSize;
    return BoundingBox<int> { colMin, rowMin, colMax, rowMax };
}

/* Resize the grid map given the bounding box (index range) */
template <typename T>
void GridMap<T>::Resize(const BoundingBox<int>& boundingBox)
{
    Assert(boundingBox.mMin.mX < boundingBox.mMax.mX);
    Assert(boundingBox.mMin.mY < boundingBox.mMax.mY);

    /* Convert the grid index to the block index */
    const Point2D<int> blockIdxMin =
        this->IndexToBlock(boundingBox.mMin.mY, boundingBox.mMin.mX);
    const Point2D<int> blockIdxMax =
        this->IndexToBlock(boundingBox.mMax.mY + this->mBlockSize - 1,
                           boundingBox.mMax.mX + this->mBlockSize - 1);

    /* Compute the number of the block rows and columns */
    const int blockRows = blockIdxMax.mY - blockIdxMin.mY;
    const int blockCols = blockIdxMax.mX - blockIdxMin.mX;
    const int numOfBlocks = blockRows * blockCols;

    /* Reallocate the blocks */
    auto oldBlocks = std::move(this->mBlocks);
    this->mBlocks.reset(new T[numOfBlocks]);
    Assert(this->mBlocks != nullptr);

    /* Move the blocks */
    BoundingBox<int> croppedBox {
        0, 0, this->mBlockCols, this->mBlockRows };
    croppedBox.Intersect(blockIdxMin.mX, blockIdxMin.mY,
                         blockIdxMax.mX, blockIdxMax.mY);

    for (int row = croppedBox.mMin.mY; row < croppedBox.mMax.mY; ++row) {
        for (int col = croppedBox.mMin.mX; col < croppedBox.mMax.mX; ++col) {
            const int blockIdx = (row - blockIdxMin.mY) * blockCols +
                                 (col - blockIdxMin.mX);
            const int oldBlockIdx = row * this->mBlockCols + col;
            this->mBlocks[blockIdx] = std::move(oldBlocks[oldBlockIdx]);
        }
    }

    /* Update the grid map parameters */
    this->mBlockRows = blockRows;
    this->mBlockCols = blockCols;

    const int rowMin = blockIdxMin.mY << this->mLog2BlockSize;
    const int colMin = blockIdxMin.mX << this->mLog2BlockSize;
    const int rows = blockRows << this->mLog2BlockSize;
    const int cols = blockCols << this->mLog2BlockSize;
    this->mGeometry.Resize(rowMin, colMin, rows, cols);
}

/* Resize the grid map given the bounding box */
template <typename T>
void GridMap<T>::Resize(const BoundingBox<double>& boundingBox)
{
    Assert(boundingBox.mMin.mX < boundingBox.mMax.mX);
    Assert(boundingBox.mMin.mY < boundingBox.mMax.mY);

    /* Avoid the rounding errors */
    const Point2D<int> idxMin = this->PositionToIndex(
        boundingBox.mMin.mX - this->mGeometry.mResolution,
        boundingBox.mMin.mY - this->mGeometry.mResolution);
    const Point2D<int> idxMax = this->PositionToIndex(
        boundingBox.mMax.mX + this->mGeometry.mResolution,
        boundingBox.mMax.mY + this->mGeometry.mResolution);

    /* Create the bounding box using the above indices */
    const BoundingBox<int> idxBoundingBox {
        idxMin.mX, idxMin.mY, idxMax.mX + 1, idxMax.mY + 1 };
    /* Resize the grid map */
    this->Resize(idxBoundingBox);
}

/* Expand the grid map given the bounding box (index range) */
template <typename T>
void GridMap<T>::Expand(const BoundingBox<int>& boundingBox)
{
    Assert(boundingBox.mMin.mX < boundingBox.mMax.mX);
    Assert(boundingBox.mMin.mY < boundingBox.mMax.mY);

    if (this->IsIndexInside(boundingBox.mMin.mY, boundingBox.mMin.mX) &&
        this->IsIndexInside(boundingBox.mMax.mY - 1, boundingBox.mMax.mX - 1))
        return;

    BoundingBox<int> expandedBox {
        0, 0, this->mGeometry.mCols, this->mGeometry.mRows };
    expandedBox.Expand(boundingBox.mMin.mX, boundingBox.mMin.mY,
                       boundingBox.mMax.mX, boundingBox.mMax.mY);

    /* Resize the grid map */
    this->Resize(expandedBox);
}

/* Expand the grid map given the bounding box */
template <typename T>
void GridMap<T>::Expand(const BoundingBox<double>& boundingBox)
{
    Assert(boundingBox.mMin.mX < boundingBox.mMax.mX);
    Assert(boundingBox.mMin.mY < boundingBox.mMax.mY);

    /* Avoid the rounding errors */
    const Point2D<int> idxMin = this->PositionToIndex(
        boundingBox.mMin.mX - this->mGeometry.mResolution,
        boundingBox.mMin.mY - this->mGeometry.mResolution);
    const Point2D<int> idxMax = this->PositionToIndex(
        boundingBox.mMax.mX + this->mGeometry.mResolution,
        boundingBox.mMax.mY + this->mGeometry.mResolution);

    /* Create the bounding box using the above indices */
    const BoundingBox<int> idxBoundingBox {
        idxMin.mX, idxMin.mY, idxMax.mX + 1, idxMax.mY + 1 };
    /* Expand the grid map */
    this->Expand(idxBoundingBox);
}

/* Crop the grid map (remove the unused blocks and shrink to fit) */
template <typename T>
void GridMap<T>::Crop()
{
    this->Resize(this->CroppedBoundingBox());
}

/* Compute the bounding box of the cropped grid map */
template <typename T>
BoundingBox<int> GridMap<T>::CroppedBoundingBox() const
{
    BoundingBox<int> croppedBox = this->CroppedBoundingBoxInBlocks();

    /* Convert the number of blocks to the number of grid cells */
    croppedBox.mMin.mX <<= this->mLog2BlockSize;
    croppedBox.mMin.mY <<= this->mLog2BlockSize;
    croppedBox.mMax.mX <<= this->mLog2BlockSize;
    croppedBox.mMax.mY <<= this->mLog2BlockSize;

    return croppedBox;
}

/* Compute the bounding box of the cropped grid map in blocks */
template <typename T>
BoundingBox<int> GridMap<T>::CroppedBoundingBoxInBlocks() const
{
    const T* block = this->Block();
    BoundingBox<int> croppedBox = BoundingBox<int>::Zero;

    /* Compute the bounding box of the allocated blocks */
    for (int row = 0; row < this->mBlockRows; ++row)
        for (int col = 0; col < this->mBlockCols; ++col, ++block)
            if (block->IsAllocated())
                croppedBox.Expand(col, row);

    croppedBox.mMax.mX += 1;
    croppedBox.mMax.mY += 1;

    return croppedBox;
}

/* Inspect the memory usage in bytes */
template <typename T>
std::uint64_t GridMap<T>::InspectMemoryUsage() const
{
    std::uint64_t memoryUsage = 0;

    memoryUsage += sizeof(this->mLog2BlockSize) +
                   sizeof(this->mBlockSize) +
                   sizeof(this->mBlockRows) +
                   sizeof(this->mBlockCols) +
                   sizeof(this->mBlocks) +
                   sizeof(this->mGeometry);

    const T* block = this->Block();

    /* Compute the memory usage for blocks */
    for (int row = 0; row < this->mBlockRows; ++row)
        for (int col = 0; col < this->mBlockCols; ++col, ++block)
            memoryUsage += block->InspectMemoryUsage();

    return memoryUsage;
}

} /* namespace GridMapNew */
} /* namespace MyGMapping */
