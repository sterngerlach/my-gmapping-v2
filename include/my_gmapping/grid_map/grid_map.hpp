
/* grid_map.hpp */

#ifndef MY_GMAPPING_GRID_MAP_GRID_MAP_HPP
#define MY_GMAPPING_GRID_MAP_GRID_MAP_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <memory>

#include "my_gmapping/grid_map/grid_map_base.hpp"
#include "my_gmapping/grid_map/grid_map_cell.hpp"
#include "my_gmapping/grid_map/grid_map_patch.hpp"
#include "my_gmapping/point.hpp"

namespace MyGMapping {

template <typename T>
class GridMap final : public GridMapBase<
    typename T::ValueType, typename T::StorageType>
{
public:
    /* Type definitions */
    using CellType = T;
    using PatchType = Patch<T>;
    using ValueType = typename T::ValueType;
    using StorageType = typename T::StorageType;
    using ObservationType = typename T::ObservationType;

    /* Constructor with number of cells */
    GridMap(double cellSize,
            int initialNumCellsX,
            int initialNumCellsY,
            const Point2D<double>& centerPos = Point2D<double>(0.0, 0.0),
            int patchSize = (1 << 4));
    
    /* Constructor with map size */
    GridMap(double cellSize,
            double mapSizeX,
            double mapSizeY,
            const Point2D<double>& centerPos = Point2D<double>(0.0, 0.0),
            int patchSize = (1 << 4));
    
    /* Constructor with positions */
    GridMap(double cellSize,
            double minX,
            double maxX,
            double minY,
            double maxY,
            int patchSize = (1 << 4));

    /* Constructor with internal parameters */
    GridMap(const double cellSize,
            const int patchSize,
            const int numOfPatchesX,
            const int numOfPatchesY,
            const double minX,
            const double minY);

    /* Destructor */
    ~GridMap() = default;

    /* Copy constructor */
    GridMap(const GridMap& other);
    /* Copy assignment operator */
    GridMap& operator=(const GridMap& other);
    /* Move constructor */
    GridMap(GridMap&& other) noexcept;
    /* Move assignment operator */
    GridMap& operator=(GridMap&& other) noexcept;

    /* Create a new grid map that has the same size as 'gridMap' */
    template <typename U>
    static GridMap CreateSameSizeMap(const GridMap<U>& gridMap);

    /* Resize the map */
    void Resize(double minX, double maxX,
                double minY, double maxY);
    /* Expand the map if necessary */
    void Expand(double minX, double maxX,
                double minY, double maxY,
                double enlargeStep = 0.1);

    /* Reset all cells inside the grid map */
    void Reset();

    /* Get the unknown occupancy probability value */
    ValueType UnknownValue() const override
    { return CellType::Unknown; }
    /* Get the internal value to represent unknown occupancy probability */
    StorageType UnknownRawValue() const override
    { return CellType::UnknownRaw; }
    /* Get the minimum internal value */
    StorageType RawValueMin() const override
    { return CellType::ValueMin; }
    /* Get the maximum internal value */
    StorageType RawValueMax() const override
    { return CellType::ValueMax; }
    /* Get the internal value range */
    StorageType RawValueRange() const override
    { return CellType::ValueRange; }

    /* Check if the index is inside the map */
    bool IsInside(int cellX, int cellY) const override;
    /* Check if the index is inside the map */
    bool IsInside(const Point2D<int>& cellIdx) const override
    { return this->IsInside(cellIdx.mX, cellIdx.mY); }

    /* Check if the point is inside the map */
    bool IsInside(double mapX, double mapY) const override
    { return this->IsInside(this->MapCoordinateToCellIndex(mapX, mapY)); }
    /* Check if the point is inside the map */
    bool IsInside(const Point2D<double>& mapPos) const override
    { return this->IsInside(mapPos.mX, mapPos.mY); }

    /* Check if the cell is allocated */
    bool IsAllocated(int cellX, int cellY) const override
    { return this->PatchIsAllocated(
        this->CellIndexToPatchIndex(cellX, cellY)); }
    /* Check if the cell is allocated */
    bool IsAllocated(const Point2D<int>& cellIdx) const override
    { return this->IsAllocated(cellIdx.mX, cellIdx.mY); }

    /* Convert cell index to map coordinate
     * Calculate the minimum position of the specified cell (x, y) */
    Point2D<double> CellIndexToMapCoordinate(int x, int y) const override;
    /* Convert cell index to map coordinate
     * Calculate the minimum position of the specified cell (x, y) */
    Point2D<double> CellIndexToMapCoordinate(
        const Point2D<int>& cellIdx) const override
    { return this->CellIndexToMapCoordinate(cellIdx.mX, cellIdx.mY); }

    /* Convert map coordinate to cell index */
    Point2D<int> MapCoordinateToCellIndex(
        const double x, const double y) const override;
    /* Convert map coordinate to cell index */
    Point2D<int> MapCoordinateToCellIndex(
        const Point2D<double>& mapPos) const override
    { return this->MapCoordinateToCellIndex(mapPos.mX, mapPos.mY); }

    /* Convert map coordinate to floating-point cell index */
    Point2D<double> MapCoordinateToCellIndexFloat(
        const double x, const double y) const override;
    /* Convert map coordinate to floating-point cell index */
    Point2D<double> MapCoordinateToCellIndexFloat(
        const Point2D<double>& mapPos) const override
    { return this->MapCoordinateToCellIndexFloat(mapPos.mX, mapPos.mY); }

    /* Get the cell of the specified index */
    CellType& At(int x, int y);
    inline CellType& At(const Point2D<int>& cellIdx)
    { return this->At(cellIdx.mX, cellIdx.mY); }

    /* Get the cell of the specified index */
    const CellType& At(int x, int y) const;
    inline const CellType& At(const Point2D<int>& cellIdx) const
    { return this->At(cellIdx.mX, cellIdx.mY); }

    /* Get the cell of the specified point */
    inline CellType& At(double mapX, double mapY)
    { return this->At(this->MapCoordinateToCellIndex(mapX, mapY)); }
    inline CellType& At(const Point2D<double>& mapPos)
    { return this->At(mapPos.mX, mapPos.mY); }

    /* Get the cell of the specified point */
    inline const CellType& At(double mapX, double mapY) const
    { return this->At(this->MapCoordinateToCellIndex(mapX, mapY)); }
    inline const CellType& At(const Point2D<double>& mapPos) const
    { return this->At(mapPos.mX, mapPos.mY); }

    /* Get the value of the specified index */
    ValueType Value(int x, int y) const override;
    ValueType Value(const Point2D<int>& cellIdx) const override
    { return this->Value(cellIdx.mX, cellIdx.mY); }

    /* Get the value of the specified index
     * The default value is returned if the specified grid cell is
     * out of bounds or is not yet allocated */
    ValueType Value(int x, int y, const ValueType defaultVal) const override;
    ValueType Value(const Point2D<int>& cellIdx,
                    const ValueType defaultVal) const override
    { return this->Value(cellIdx.mX, cellIdx.mY, defaultVal); }

    /* Get the internal value of the specified grid cell */
    StorageType RawValue(const int idxX, const int idxY) const override;
    /* Get the internal value of the specified grid cell */
    StorageType RawValue(const Point2D<int>& cellIdx) const override
    { return this->RawValue(cellIdx.mX, cellIdx.mY); }

    /* Get the internal value of the specified grid cell
     * The default value is returned if the specified grid cell is
     * out of bounds or is not yet allocated */
    StorageType RawValue(const int idxX, const int idxY,
                         const StorageType defaultVal) const override;
    /* Get the internal value of the specified grid cell */
    StorageType RawValue(const Point2D<int>& cellIdx,
                         const StorageType defaultVal) const override
    { return this->RawValue(cellIdx.mX, cellIdx.mY, defaultVal); }

    /* Set the occupancy probability value of the specified grid cell */
    void SetValue(const int idxX, const int idxY, const ValueType probValue);
    /* Set the occupancy probability value of the specified grid cell */
    void SetValue(const Point2D<int>& cellIdx,
                  const ValueType probValue)
    { this->SetValue(cellIdx.mX, cellIdx.mY, probValue); }

    /* Set the internal value of the specified grid cell */
    void SetRawValue(const int idxX, const int idxY,
                     const StorageType rawValue);
    /* Set the internal value of the specified grid cell */
    void SetRawValue(const Point2D<int>& cellIdx,
                     const ValueType rawValue)
    { this->SetRawValue(cellIdx.mX, cellIdx.mY, rawValue); }

    /* Update the value of the specified cell */
    void Update(int x, int y, const ObservationType& updateVal);
    inline void Update(const Point2D<int>& cellIdx,
                       const ObservationType& updateVal)
    { this->Update(cellIdx.mX, cellIdx.mY, updateVal); }

    /* Update the value of the specified point */
    void Update(double mapX, double mapY, const ObservationType& updateVal)
    { this->Update(this->MapCoordinateToCellIndex(mapX, mapY), updateVal); }
    inline void Update(const Point2D<double>& mapPos,
                       const ObservationType& updateVal)
    { this->Update(mapPos.mX, mapPos.mY, updateVal); }

    /* Calculate the distance between two cell indices */
    double Distance(int x0, int y0, int x1, int y1) const override;
    double Distance(const Point2D<int>& cellIdx0,
                    const Point2D<int>& cellIdx1) const override
    { return this->Distance(cellIdx0.mX, cellIdx0.mY,
                            cellIdx1.mX, cellIdx1.mY); }
    
    /* Calculate the squared distance between two cell indices */
    double SquaredDistance(int x0, int y0, int x1, int y1) const override;
    double SquaredDistance(const Point2D<int>& cellIdx0,
                           const Point2D<int>& cellIdx1) const override
    { return this->SquaredDistance(cellIdx0.mX, cellIdx0.mY,
                                   cellIdx1.mX, cellIdx1.mY); }
    
    /* Convert cell index to patch index */
    Point2D<int> CellIndexToPatchIndex(int x, int y) const;
    /* Convert cell index to patch index */
    inline Point2D<int> CellIndexToPatchIndex(const Point2D<int>& cellIdx) const
    { return this->CellIndexToPatchIndex(cellIdx.mX, cellIdx.mY); }

    /* Convert the cell index to the offset inside the patch */
    Point2D<int> CellIndexToPatchOffset(int x, int y) const;
    /* Convert the cell index to the offset inside the patch */
    inline Point2D<int> CellIndexToPatchOffset(
        const Point2D<int>& cellIdx) const
    { return this->CellIndexToPatchOffset(cellIdx.mX, cellIdx.mY); }

    /* Convert patch index to cell index range */
    void PatchIndexToCellIndexRange(int x, int y,
                                    int& cellMinX, int& cellMinY,
                                    int& cellMaxX, int& cellMaxY) const;
    void PatchIndexToCellIndexRange(const Point2D<int>& patchIdx,
                                    Point2D<int>& cellMinIdx,
                                    Point2D<int>& cellMaxIdx) const
    { this->PatchIndexToCellIndexRange(patchIdx.mX, patchIdx.mY,
                                       cellMinIdx.mX, cellMinIdx.mY,
                                       cellMaxIdx.mX, cellMaxIdx.mY); }

    /* Check if the patch is inside the map */
    bool PatchIsInside(int patchIdxX, int patchIdxY) const;
    /* Check if the patch is inside the map */
    bool PatchIsInside(const Point2D<int>& patchIdx) const
    { return this->PatchIsInside(patchIdx.mX, patchIdx.mY); }

    /* Check if the patch is allocated on the heap */
    bool PatchIsAllocated(int patchIdxX, int patchIdxY) const;
    /* Check if the patch is allocated on the heap */
    bool PatchIsAllocated(const Point2D<int>& patchIdx) const
    { return this->PatchIsAllocated(patchIdx.mX, patchIdx.mY); }

    /* Get the patch of the specified index */
    inline Patch<T>& PatchAt(int patchIdxX, int patchIdxY);
    inline Patch<T>& PatchAt(const Point2D<int>& patchIdx)
    { return this->PatchAt(patchIdx.mX, patchIdx.mY); }

    /* Get the patch of the specified index */
    inline const Patch<T>& PatchAt(int patchIdxX, int patchIdxY) const;
    inline const Patch<T>& PatchAt(const Point2D<int>& patchIdx) const
    { return this->PatchAt(patchIdx.mX, patchIdx.mY); }

    /* Get the pointer to the specified patch */
    const Patch<T>* PatchPtrAt(int patchIdxX, int patchIdxY) const;
    /* Get the pointer to the specified patch */
    const Patch<T>* PatchPtrAt(const Point2D<int>& patchIdx) const
    { return this->PatchPtrAt(patchIdx.mX, patchIdx.mY); }

    /* Compute the actual grid map size */
    void ComputeActualMapSize(Point2D<int>& patchIdxMin,
                              Point2D<int>& patchIdxMax,
                              Point2D<int>& gridCellIdxMin,
                              Point2D<int>& gridCellIdxMax,
                              Point2D<int>& mapSizeInPatches,
                              Point2D<int>& mapSizeInGridCells) const;

    /* Get the cell size */
    inline double CellSize() const override { return this->mCellSize; }

    /* Get the number of the patches */
    inline int NumPatchX() const { return this->mNumPatchX; }
    inline int NumPatchY() const { return this->mNumPatchY; }

    /* Get the number of the cells */
    inline int NumCellsX() const override { return this->mNumCellsX; }
    inline int NumCellsY() const override { return this->mNumCellsY; }
    
    /* Get the size of the map */
    inline double MapSizeX() const override { return this->mMapSizeX; }
    inline double MapSizeY() const override { return this->mMapSizeY; }

    /* Get the size of the patch */
    inline int PatchSize() const { return this->mPatchSize; }
    
    /* Get the minimum position of the grid map (in world coordinate) */
    inline const Point2D<double>& MinPos() const override
    { return this->mMinPos; }
    /* Set the minimum position of the grid map (in world coordinate) */
    inline void SetMinPos(const double minX, const double minY)
    { this->mMinPos.mX = minX; this->mMinPos.mY = minY; }
    /* Set the minimum position of the grid map (in world coordinate) */
    inline void SetMinPos(const Point2D<double>& minPos)
    { this->mMinPos = minPos; }

private:
    /* Map resolution (grid cell size in meters) */
    double                      mCellSize;
    /* Number of the patches (horizontal) */
    int                         mNumPatchX;
    /* Number of the patches (vertical) */
    int                         mNumPatchY;
    /* Number of the grid cells (horizontal) */
    int                         mNumCellsX;
    /* Number of the grid cells (vertical) */
    int                         mNumCellsY;
    /* Size of the map in meters (horizontal) */
    double                      mMapSizeX;
    /* Size of the map in meters (vertical) */
    double                      mMapSizeY;
    /* Size of the patch (in the number of grid cells) */
    int                         mPatchSize;
    /* Minimum position of the grid map (in world coordinate)
     * corresponding to the top-left of the origin grid cell (0, 0) */
    Point2D<double>             mMinPos;
    /* Patches */
    std::shared_ptr<Patch<T>[]> mPatches;
};

/* Constructor with number of cells */
template <typename T>
GridMap<T>::GridMap(double cellSize,
                    int initialNumCellsX,
                    int initialNumCellsY,
                    const Point2D<double>& centerPos,
                    int patchSize)
{
    assert(cellSize > 0.0);
    assert(initialNumCellsX >= 0);
    assert(initialNumCellsY >= 0);
    assert(patchSize > 0);

    initialNumCellsX = std::max(1, initialNumCellsX);
    initialNumCellsY = std::max(1, initialNumCellsY);

    this->mCellSize = cellSize;

    this->mNumPatchX = static_cast<int>(std::ceil(
        static_cast<double>(initialNumCellsX) /
        static_cast<double>(patchSize)));
    this->mNumPatchY = static_cast<int>(std::ceil(
        static_cast<double>(initialNumCellsY) /
        static_cast<double>(patchSize)));

    assert(this->mNumPatchX > 0);
    assert(this->mNumPatchY > 0);

    this->mPatchSize = patchSize;
    this->mNumCellsX = this->mNumPatchX * this->mPatchSize;
    this->mNumCellsY = this->mNumPatchY * this->mPatchSize;
    this->mMapSizeX = static_cast<double>(this->mNumCellsX) * cellSize;
    this->mMapSizeY = static_cast<double>(this->mNumCellsY) * cellSize;

    const double idxOffsetX = (this->mNumCellsX % 2 == 0) ?
        (this->mNumCellsX / 2) : (this->mNumCellsX / 2 + 0.5);
    const double idxOffsetY = (this->mNumCellsY % 2 == 0) ?
        (this->mNumCellsY / 2) : (this->mNumCellsY / 2 + 0.5);
    this->mMinPos.mX = centerPos.mX - idxOffsetX * this->mCellSize;
    this->mMinPos.mY = centerPos.mY - idxOffsetY * this->mCellSize;

    this->mPatches.reset(new Patch<T>[this->mNumPatchX * this->mNumPatchY]);
    assert(this->mPatches != nullptr);
}

/* Constructor with map size */
template <typename T>
GridMap<T>::GridMap(double cellSize,
                    double mapSizeX,
                    double mapSizeY,
                    const Point2D<double>& centerPos,
                    int patchSize)
{
    assert(cellSize > 0.0);
    assert(mapSizeX >= 0.0);
    assert(mapSizeY >= 0.0);
    assert(patchSize > 0);

    const double patchSizeInMeters =
        static_cast<double>(patchSize) * cellSize;

    mapSizeX = std::max(patchSizeInMeters, mapSizeX);
    mapSizeY = std::max(patchSizeInMeters, mapSizeY);

    this->mCellSize = cellSize;

    this->mNumPatchX = static_cast<int>(std::ceil(
        mapSizeX / patchSizeInMeters));
    this->mNumPatchY = static_cast<int>(std::ceil(
        mapSizeY / patchSizeInMeters));

    assert(this->mNumPatchX > 0);
    assert(this->mNumPatchY > 0);

    this->mPatchSize = patchSize;
    this->mNumCellsX = this->mNumPatchX * this->mPatchSize;
    this->mNumCellsY = this->mNumPatchY * this->mPatchSize;
    this->mMapSizeX = static_cast<double>(this->mNumCellsX) * cellSize;
    this->mMapSizeY = static_cast<double>(this->mNumCellsY) * cellSize;

    const double idxOffsetX = (this->mNumCellsX % 2 == 0) ?
        (this->mNumCellsX / 2) : (this->mNumCellsX / 2 + 0.5);
    const double idxOffsetY = (this->mNumCellsY % 2 == 0) ?
        (this->mNumCellsY / 2) : (this->mNumCellsY / 2 + 0.5);
    this->mMinPos.mX = centerPos.mX - idxOffsetX * this->mCellSize;
    this->mMinPos.mY = centerPos.mY - idxOffsetY * this->mCellSize;

    this->mPatches.reset(new Patch<T>[this->mNumPatchX * this->mNumPatchY]);
    assert(this->mPatches != nullptr);
}

/* Constructor with positions */
template <typename T>
GridMap<T>::GridMap(double cellSize,
                    double minX,
                    double maxX,
                    double minY,
                    double maxY,
                    int patchSize)
{
    assert(cellSize > 0.0);
    assert(minX <= maxX);
    assert(minY <= maxY);
    assert(patchSize > 0);

    const double patchSizeInMeters =
        static_cast<double>(patchSize) * cellSize;
    double mapSizeX = std::max(patchSizeInMeters, maxX - minX);
    double mapSizeY = std::max(patchSizeInMeters, maxY - minY);

    this->mCellSize = cellSize;

    this->mNumPatchX = static_cast<int>(std::ceil(
        mapSizeX / patchSizeInMeters));
    this->mNumPatchY = static_cast<int>(std::ceil(
        mapSizeY / patchSizeInMeters));

    assert(this->mNumPatchX > 0);
    assert(this->mNumPatchY > 0);

    this->mPatchSize = patchSize;
    this->mNumCellsX = this->mNumPatchX * this->mPatchSize;
    this->mNumCellsY = this->mNumPatchY * this->mPatchSize;
    this->mMapSizeX = static_cast<double>(this->mNumCellsX) * cellSize;
    this->mMapSizeY = static_cast<double>(this->mNumCellsY) * cellSize;

    this->mMinPos.mX = minX;
    this->mMinPos.mY = minY;

    this->mPatches.reset(new Patch<T>[this->mNumPatchX * this->mNumPatchY]);
    assert(this->mPatches != nullptr);
}

/* Constructor with internal parameters */
template <typename T>
GridMap<T>::GridMap(const double cellSize,
                    const int patchSize,
                    const int numOfPatchesX,
                    const int numOfPatchesY,
                    const double minX,
                    const double minY)
{
    /* Input validity checks */
    assert(cellSize > 0.0);
    assert(patchSize > 0);
    assert(numOfPatchesX > 0);
    assert(numOfPatchesY > 0);

    /* Set the map resolution */
    this->mCellSize = cellSize;
    /* Set the patch size */
    this->mPatchSize = patchSize;

    /* Set the number of patches */
    this->mNumPatchX = numOfPatchesX;
    this->mNumPatchY = numOfPatchesY;
    /* Set the number of grid cells */
    this->mNumCellsX = this->mNumPatchX * this->mPatchSize;
    this->mNumCellsY = this->mNumPatchY * this->mPatchSize;
    /* Set the map size in meters */
    this->mMapSizeX = this->mNumCellsX * this->mCellSize;
    this->mMapSizeY = this->mNumCellsY * this->mCellSize;
    /* Set the minimum position */
    this->mMinPos.mX = minX;
    this->mMinPos.mY = minY;

    /* Allocate the memory for patches */
    this->mPatches.reset(new Patch<T>[this->mNumPatchX * this->mNumPatchY]);
    assert(this->mPatches != nullptr);
}

/* Copy constructor */
template <typename T>
GridMap<T>::GridMap(const GridMap<T>& other) :
    GridMapBase<typename T::ValueType,
                typename T::StorageType>(other),
    mCellSize(other.mCellSize),
    mNumPatchX(other.mNumPatchX),
    mNumPatchY(other.mNumPatchY),
    mNumCellsX(other.mNumCellsX),
    mNumCellsY(other.mNumCellsY),
    mMapSizeX(other.mMapSizeX),
    mMapSizeY(other.mMapSizeY),
    mPatchSize(other.mPatchSize),
    mMinPos(other.mMinPos),
    mPatches(nullptr)
{
    /* Allocate patches then copy the values */
    assert(this->mNumPatchX > 0);
    assert(this->mNumPatchY > 0);

    this->mPatches.reset(new Patch<T>[this->mNumPatchX * this->mNumPatchY]);
    assert(this->mPatches != nullptr);

    /* Copy the patches */
    std::copy_n(other.mPatches.get(), this->mNumPatchX * this->mNumPatchY,
                this->mPatches.get());
}

/* Copy assignment operator */
template <typename T>
GridMap<T>& GridMap<T>::operator=(const GridMap<T>& other)
{
    if (this == &other)
        return *this;
    
    /* Copy the grid map parameters */
    this->mCellSize = other.mCellSize;
    this->mNumPatchX = other.mNumPatchX;
    this->mNumPatchY = other.mNumPatchY;
    this->mNumCellsX = other.mNumCellsX;
    this->mNumCellsY = other.mNumCellsY;
    this->mMapSizeX = other.mMapSizeX;
    this->mMapSizeY = other.mMapSizeY;
    this->mPatchSize = other.mPatchSize;
    this->mMinPos = other.mMinPos;

    /* Allocate patches then copy the values */
    assert(this->mNumPatchX > 0);
    assert(this->mNumPatchY > 0);

    this->mPatches.reset(new Patch<T>[this->mNumPatchX * this->mNumPatchY]);
    assert(this->mPatches != nullptr);

    /* Copy the patches */
    std::copy_n(other.mPatches.get(), this->mNumPatchX * this->mNumPatchY,
                this->mPatches.get());
    
    return *this;
}

/* Move constructor */
template <typename T>
GridMap<T>::GridMap(GridMap<T>&& other) noexcept :
    mCellSize(other.mCellSize),
    mNumPatchX(other.mNumPatchX),
    mNumPatchY(other.mNumPatchY),
    mNumCellsX(other.mNumCellsX),
    mNumCellsY(other.mNumCellsY),
    mMapSizeX(other.mMapSizeX),
    mMapSizeY(other.mMapSizeY),
    mPatchSize(other.mPatchSize),
    mMinPos(other.mMinPos),
    mPatches(std::move(other.mPatches))
{
}

/* Move assignment operator */
template <typename T>
GridMap<T>& GridMap<T>::operator=(GridMap<T>&& other) noexcept
{
    if (this == &other)
        return *this;
    
    /* Copy the grid map parameters */
    this->mCellSize = other.mCellSize;
    this->mNumPatchX = other.mNumPatchX;
    this->mNumPatchY = other.mNumPatchY;
    this->mNumCellsX = other.mNumCellsX;
    this->mNumCellsY = other.mNumCellsY;
    this->mMapSizeX = other.mMapSizeX;
    this->mMapSizeY = other.mMapSizeY;
    this->mPatchSize = other.mPatchSize;
    this->mMinPos = other.mMinPos;

    /* Move and assign the patches */
    this->mPatches = std::move(other.mPatches);

    return *this;
}

/* Create a new grid map that has the same size as 'gridMap' */
template <typename T>
template <typename U>
GridMap<T> GridMap<T>::CreateSameSizeMap(const GridMap<U>& gridMap)
{
    /* Construct a new grid map with internal parameters */
    return GridMap<T> { gridMap.CellSize(), gridMap.PatchSize(),
                        gridMap.NumPatchX(), gridMap.NumPatchY(),
                        gridMap.MinPos().mX, gridMap.MinPos().mY };
}

/* Resize the map */
template <typename T>
void GridMap<T>::Resize(double minX, double maxX,
                        double minY, double maxY)
{
    assert(minX <= maxX);
    assert(minY <= maxY);

    /* Calculate the cell index */
    const Point2D<int> cellMinIdx = this->MapCoordinateToCellIndex(minX, minY);
    const Point2D<int> cellMaxIdx = this->MapCoordinateToCellIndex(maxX, maxY);

    /* Calculate the patch index (allow negative index) */
    const Point2D<int> patchMinIdx = this->CellIndexToPatchIndex(cellMinIdx);
    const Point2D<int> patchMaxIdx = this->CellIndexToPatchIndex(cellMaxIdx);

    /* Calculate the number of patches */
    const int newNumPatchX = patchMaxIdx.mX - patchMinIdx.mX + 1;
    const int newNumPatchY = patchMaxIdx.mY - patchMinIdx.mY + 1;

    /* Reallocate patches */
    std::shared_ptr<Patch<T>[]> oldPatches = std::move(this->mPatches);
    this->mPatches.reset(new Patch<T>[newNumPatchX * newNumPatchY]);
    assert(this->mPatches != nullptr);

    /* Copy patches */
    const int x0 = std::max(0, patchMinIdx.mX);
    const int y0 = std::max(0, patchMinIdx.mY);
    const int x1 = std::min(this->mNumPatchX, patchMaxIdx.mX + 1);
    const int y1 = std::min(this->mNumPatchY, patchMaxIdx.mY + 1);

    for (int y = y0; y < y1; ++y) {
        for (int x = x0; x < x1; ++x) {
            const int patchIdx = (y - patchMinIdx.mY) * newNumPatchX +
                                 (x - patchMinIdx.mX);
            const int oldPatchIdx = y * this->mNumPatchX + x;
            this->mPatches[patchIdx] = std::move(oldPatches[oldPatchIdx]);
        }
    }

    /* Update map parameters */
    this->mNumPatchX = newNumPatchX;
    this->mNumPatchY = newNumPatchY;
    this->mNumCellsX = this->mNumPatchX * this->mPatchSize;
    this->mNumCellsY = this->mNumPatchY * this->mPatchSize;
    this->mMapSizeX = this->mNumCellsX * this->mCellSize;
    this->mMapSizeY = this->mNumCellsY * this->mCellSize;

    this->mMinPos.mX += (patchMinIdx.mX * this->mPatchSize) * this->mCellSize;
    this->mMinPos.mY += (patchMinIdx.mY * this->mPatchSize) * this->mCellSize;
}

/* Expand the map if necessary */
template <typename T>
void GridMap<T>::Expand(double minX, double maxX,
                        double minY, double maxY,
                        double enlargeStep)
{
    assert(minX <= maxX);
    assert(minY <= maxY);

    if (this->IsInside(minX, minY) && this->IsInside(maxX, maxY))
        return;

    /* Retrieve the minimum position of the grid map */
    Point2D<double> minPos = this->mMinPos;
    /* Calculate the maximum position of the grid map */
    Point2D<double> maxPos = this->CellIndexToMapCoordinate(
        this->mNumCellsX, this->mNumCellsY);

    minPos.mX = (minX < minPos.mX) ? minX - enlargeStep : minPos.mX;
    minPos.mY = (minY < minPos.mY) ? minY - enlargeStep : minPos.mY;
    maxPos.mX = (maxX > maxPos.mX) ? maxX + enlargeStep : maxPos.mX;
    maxPos.mY = (maxY > maxPos.mY) ? maxY + enlargeStep : maxPos.mY;

    /* Expand the map if necessary */
    this->Resize(minPos.mX, maxPos.mX, minPos.mY, maxPos.mY);
}

/* Reset all cells inside the grid map */
template <typename T>
void GridMap<T>::Reset()
{
    assert(this->mNumPatchX > 0);
    assert(this->mNumPatchY > 0);
    assert(this->mPatches != nullptr);

    /* Reset all patches and cells */
    for (int patchIdxY = 0; patchIdxY < this->mNumPatchY; ++patchIdxY)
        for (int patchIdxX = 0; patchIdxX < this->mNumPatchX; ++patchIdxX)
            this->PatchAt(patchIdxX, patchIdxY).Reset();
}

/* Check if the index is inside the map */
template <typename T>
bool GridMap<T>::IsInside(int cellX, int cellY) const
{
    return (cellX >= 0 && cellX < this->mNumCellsX) &&
           (cellY >= 0 && cellY < this->mNumCellsY);
}

/* Convert cell index to map coordinate
 * Calculate the minimum position of the specified cell (x, y) */
template <typename T>
Point2D<double> GridMap<T>::CellIndexToMapCoordinate(int x, int y) const
{
    /* Calculate the minimum position of the specified cell (x, y)
     * this->mMinPos represents the minimum position of the map */
    const double mapX = this->mMinPos.mX + this->mCellSize * x;
    const double mapY = this->mMinPos.mY + this->mCellSize * y;

    return Point2D<double> { mapX, mapY };
}

/* Convert map coordinate to cell index */
template <typename T>
Point2D<int> GridMap<T>::MapCoordinateToCellIndex(
    const double x, const double y) const
{
    const int cellX = static_cast<int>(std::floor(
        (x - this->mMinPos.mX) / this->mCellSize));
    const int cellY = static_cast<int>(std::floor(
        (y - this->mMinPos.mY) / this->mCellSize));

    return Point2D<int> { cellX, cellY };
}

/* Convert map coordinate to floating-point cell index */
template <typename T>
Point2D<double> GridMap<T>::MapCoordinateToCellIndexFloat(
    const double x, const double y) const
{
    const double cellX = (x - this->mMinPos.mX) / this->mCellSize;
    const double cellY = (y - this->mMinPos.mY) / this->mCellSize;

    return Point2D<double> { cellX, cellY };
}

/* Get the cell of the specified index */
template <typename T>
typename GridMap<T>::CellType& GridMap<T>::At(int x, int y)
{
    assert(this->IsInside(x, y));

    const Point2D<int> patchIdx = this->CellIndexToPatchIndex(x, y);
    Patch<T>& patch = this->PatchAt(patchIdx);

    /* Allocate patch if necessary */
    if (!patch.IsAllocated())
        patch.Allocate(this->mPatchSize);

    return patch.At(x % this->mPatchSize,
                    y % this->mPatchSize);
}

/* Get the cell of the specified index */
template <typename T>
const typename GridMap<T>::CellType& GridMap<T>::At(int x, int y) const
{
    assert(this->IsInside(x, y));

    const Point2D<int> patchIdx = this->CellIndexToPatchIndex(x, y);
    const Patch<T>& patch = this->PatchAt(patchIdx);
    
    return patch.At(x % this->mPatchSize,
                    y % this->mPatchSize);
}

/* Get the value of the specified index */
template <typename T>
typename GridMap<T>::ValueType GridMap<T>::Value(int x, int y) const
{
    assert(this-IsInside(x, y));

    const Point2D<int> patchIdx = this->CellIndexToPatchIndex(x, y);
    const Patch<T>& patch = this->PatchAt(patchIdx);

    return patch.Value(x % this->mPatchSize,
                       y % this->mPatchSize);
}

/* Get the value of the specified index */
template <typename T>
typename GridMap<T>::ValueType GridMap<T>::Value(
    int x, int y, const ValueType defaultVal) const
{
    if (!this->IsInside(x, y))
        return defaultVal;

    const Point2D<int> patchIdx = this->CellIndexToPatchIndex(x, y);
    const Patch<T>& patch = this->PatchAt(patchIdx);

    return patch.Value(x % this->mPatchSize,
                       y % this->mPatchSize,
                       defaultVal);
}

/* Get the internal value of the specified grid cell */
template <typename T>
typename GridMap<T>::StorageType GridMap<T>::RawValue(
    const int idxX, const int idxY) const
{
    assert(this-IsInside(idxX, idxY));

    const Point2D<int> patchIdx = this->CellIndexToPatchIndex(idxX, idxY);
    const Patch<T>& patch = this->PatchAt(patchIdx);

    return patch.RawValue(idxX % this->mPatchSize,
                          idxY % this->mPatchSize);
}

/* Get the internal value of the specified grid cell */
template <typename T>
typename GridMap<T>::StorageType GridMap<T>::RawValue(
    const int idxX, const int idxY, const StorageType defaultVal) const
{
    if (!this->IsInside(idxX, idxY))
        return defaultVal;

    const Point2D<int> patchIdx = this->CellIndexToPatchIndex(idxX, idxY);
    const Patch<T>& patch = this->PatchAt(patchIdx);

    return patch.RawValue(idxX % this->mPatchSize,
                          idxY % this->mPatchSize,
                          defaultVal);
}

/* Set the occupancy probability value of the specified grid cell */
template <typename T>
void GridMap<T>::SetValue(const int idxX, const int idxY,
                          const ValueType probValue)
{
    this->At(idxX, idxY).SetValue(probValue);
}

/* Set the internal value of the specified grid cell */
template <typename T>
void GridMap<T>::SetRawValue(const int idxX, const int idxY,
                             const StorageType rawValue)
{
    this->At(idxX, idxY).SetRawValue(rawValue);
}

/* Update the value of the specified cell */
template <typename T>
void GridMap<T>::Update(int x, int y, const ObservationType& updateVal)
{
    this->At(x, y).Update(updateVal);
}

/* Calculate the distance between two cell indices */
template <typename T>
double GridMap<T>::Distance(int x0, int y0, int x1, int y1) const
{
    return this->mCellSize *
           std::hypot(static_cast<double>(x1 - x0),
                      static_cast<double>(y1 - y0));
}

/* Calculate the squared distance between two cell indices */
template <typename T>
double GridMap<T>::SquaredDistance(int x0, int y0, int x1, int y1) const
{
    const double deltaX = (x1 - x0) * this->mCellSize;
    const double deltaY = (y1 - y0) * this->mCellSize;
    return deltaX * deltaX + deltaY * deltaY;
}

/* Convert cell index to patch index */
template <typename T>
Point2D<int> GridMap<T>::CellIndexToPatchIndex(int x, int y) const
{
    /* Cell index could be negative
     * Could also be implemented using std::floor() */
    return Point2D<int> {
        (x < 0) ? (x / this->mPatchSize - 1) : (x / this->mPatchSize),
        (y < 0) ? (y / this->mPatchSize - 1) : (y / this->mPatchSize) };
}

/* Convert the cell index to the offset inside the patch */
template <typename T>
Point2D<int> GridMap<T>::CellIndexToPatchOffset(int x, int y) const
{
    return Point2D<int> { x % this->mPatchSize,
                          y % this->mPatchSize };
}

/* Convert patch index to cell index range */
template <typename T>
void GridMap<T>::PatchIndexToCellIndexRange(int x, int y,
                                            int& cellMinX, int& cellMinY,
                                            int& cellMaxX, int& cellMaxY) const
{
    cellMinX = x * this->mPatchSize;
    cellMinY = y * this->mPatchSize;
    cellMaxX = (x + 1) * this->mPatchSize;
    cellMaxY = (y + 1) * this->mPatchSize;
}

/* Check if the patch is inside the map */
template <typename T>
bool GridMap<T>::PatchIsInside(int patchIdxX, int patchIdxY) const
{
    return (patchIdxX >= 0 && patchIdxX < this->mNumPatchX) &&
           (patchIdxY >= 0 && patchIdxY < this->mNumPatchY);
}

/* Check if the patch is allocated on the heap */
template <typename T>
bool GridMap<T>::PatchIsAllocated(int patchIdxX, int patchIdxY) const
{
    assert(this->PatchIsInside(patchIdxX, patchIdxY));

    const int patchIdx = patchIdxY * this->mNumPatchX + patchIdxX;
    const Patch<T>& patch = this->mPatches[patchIdx];

    return patch.IsAllocated();
}

/* Get the patch of the specified index */
template <typename T>
Patch<T>& GridMap<T>::PatchAt(int patchIdxX, int patchIdxY)
{
    assert(this->PatchIsInside(patchIdxX, patchIdxY));
    return this->mPatches[patchIdxY * this->mNumPatchX + patchIdxX];
}

/* Get the patch of the specified index */
template <typename T>
const Patch<T>& GridMap<T>::PatchAt(int patchIdxX, int patchIdxY) const
{
    assert(this->PatchIsInside(patchIdxX, patchIdxY));
    return this->mPatches[patchIdxY * this->mNumPatchX + patchIdxX];
}

/* Get the pointer to the specified patch */
template <typename T>
const Patch<T>* GridMap<T>::PatchPtrAt(int patchIdxX, int patchIdxY) const
{
    assert(this->PatchIsInside(patchIdxX, patchIdxX));
    return this->mPatches.get() + patchIdxY * this->mNumPatchX + patchIdxX;
}

/* Compute the actual grid map size */
template <typename T>
void GridMap<T>::ComputeActualMapSize(Point2D<int>& patchIdxMin,
                                      Point2D<int>& patchIdxMax,
                                      Point2D<int>& gridCellIdxMin,
                                      Point2D<int>& gridCellIdxMax,
                                      Point2D<int>& mapSizeInPatches,
                                      Point2D<int>& mapSizeInGridCells) const
{
    /* Determine the range of the grid patch index */
    patchIdxMin.mX = std::numeric_limits<int>::max();
    patchIdxMin.mY = std::numeric_limits<int>::max();
    patchIdxMax.mX = std::numeric_limits<int>::min();
    patchIdxMax.mY = std::numeric_limits<int>::min();

    for (int y = 0; y < this->mNumPatchY; ++y) {
        for (int x = 0; x < this->mNumPatchX; ++x) {
            if (!this->PatchIsAllocated(x, y))
                continue;

            patchIdxMin.mX = std::min(patchIdxMin.mX, x);
            patchIdxMin.mY = std::min(patchIdxMin.mY, y);
            patchIdxMax.mX = std::max(patchIdxMax.mX, x);
            patchIdxMax.mY = std::max(patchIdxMax.mY, y);
        }
    }

    /* Determine the range of the grid cell index */
    Point2D<int> gridCellIdxTmp;
    this->PatchIndexToCellIndexRange(
        patchIdxMin, gridCellIdxMin, gridCellIdxTmp);
    this->PatchIndexToCellIndexRange(
        patchIdxMax, gridCellIdxTmp, gridCellIdxMax);

    /* Correct the maximum grid patch index */
    patchIdxMax.mX += 1;
    patchIdxMax.mY += 1;

    /* Compute the actual grid map size */
    mapSizeInPatches.mX = patchIdxMax.mX - patchIdxMin.mX;
    mapSizeInPatches.mY = patchIdxMax.mY - patchIdxMin.mY;
    mapSizeInGridCells.mX = mapSizeInPatches.mX * this->mPatchSize;
    mapSizeInGridCells.mY = mapSizeInPatches.mY * this->mPatchSize;

    assert(mapSizeInGridCells.mX == gridCellIdxMax.mX - gridCellIdxMin.mX);
    assert(mapSizeInGridCells.mY == gridCellIdxMax.mY - gridCellIdxMin.mY);
}

} /* namespace MyGMapping */

#endif /* MY_GMAPPING_GRID_MAP_GRID_MAP_HPP */
