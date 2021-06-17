
/* grid_map_patch.hpp */

#ifndef MY_GMAPPING_GRID_MAP_GRID_MAP_PATCH_HPP
#define MY_GMAPPING_GRID_MAP_GRID_MAP_PATCH_HPP

#include <algorithm>
#include <cassert>
#include <memory>
#include <utility>

namespace MyGMapping {

template <typename T>
class Patch final
{
public:
    /* Type definitions */
    using ValueType = typename T::ValueType;
    using StorageType = typename T::StorageType;
    using CellType = T;

    /* Default constructor */
    Patch() : mSize(-1), mData(nullptr), mIsUnique(true) { }
    /* Destructor */
    ~Patch() = default;

    /* Copy constructor */
    Patch(const Patch& other);
    /* Copy assignment operator */
    Patch& operator=(const Patch& other);

    /* Move constructor */
    Patch(Patch&& other) noexcept;
    /* Move assignment operator */
    Patch& operator=(Patch&& other) noexcept;

    /* Allocate cells */
    void Allocate(int patchSize);

    /* Check if patch is allocated */
    inline bool IsAllocated() const { return this->mData != nullptr; }

    /* Get the size of the patch */
    inline int Size() const { return this->mSize; }

    /* Get the cell of the specified index inside this patch */
    CellType& At(int x, int y);
    /* Get the cell of the specified index inside this patch */
    const CellType& At(int x, int y) const;

    /* Get the cell value of the patch */
    const ValueType Value(int x, int y) const;
    /* Get the cell value of the patch
     * Return the default value if not allocated */
    const ValueType Value(int x, int y, const ValueType defaultVal) const;

    /* Get the internal value of the specified grid cell */
    StorageType RawValue(const int x, const int y) const;
    /* Get the internal value of the specified grid cell */
    inline StorageType RawValue(const Point2D<int>& gridCellIdx) const
    { return this->RawValue(gridCellIdx.mX, gridCellIdx.mY); }

    /* Get the internal value of the specified grid cell */
    StorageType RawValue(const int x, const int y,
                         const StorageType defaultVal) const;
    /* Get the internal value of the specified grid cell */
    inline StorageType RawValue(const Point2D<int>& gridCellIdx,
                                const StorageType defaultVal) const
    { return this->RawValue(gridCellIdx.mX, gridCellIdx.mY, defaultVal); }

    /* Reset all cells inside the patch */
    void Reset();

private:
    /* Validate the patch */
    inline bool IsValid() const;
    /* Check if the index is inside the patch */
    inline bool IsInside(int cellX, int cellY) const;

private:
    /* Patch size */
    int mSize;

    /* Patch data */
    std::shared_ptr<T[]> mData;

    /* Flag to determine whether the patch is unique
     * (not being shared among several particles) */
    /* The grid cells in the patch is shared among
     * several particles to reduce memory consumption and
     * perform resampling process efficiently */
    bool mIsUnique;
};

/* Copy constructor */
template <typename T>
Patch<T>::Patch(const Patch<T>& other) :
    mSize(other.mSize),
    mData(other.mData),
    mIsUnique(false)
{
}

/* Copy assignment operator */
template <typename T>
Patch<T>& Patch<T>::operator=(const Patch<T>& other)
{
    if (this == &other)
        return *this;
    
    this->mSize = other.mSize;
    this->mData = other.mData;
    this->mIsUnique = false;

    return *this;
}

/* Move constructor */
template <typename T>
Patch<T>::Patch(Patch<T>&& other) noexcept :
    mSize(other.mSize),
    mData(std::move(other.mData)),
    mIsUnique(other.mIsUnique)
{
}

/* Move assignment operator */
template <typename T>
Patch<T>& Patch<T>::operator=(Patch<T>&& other) noexcept
{
    if (this == &other)
        return *this;
    
    this->mSize = other.mSize;
    this->mData = std::move(other.mData);
    this->mIsUnique = other.mIsUnique;

    return *this;
}

/* Allocate patch data */
template <typename T>
void Patch<T>::Allocate(int patchSize)
{    
    assert(patchSize > 0);

    this->mSize = patchSize;
    this->mData.reset(new T[patchSize * patchSize]);
    this->mIsUnique = true;

    assert(this->mData != nullptr);
}

/* Get the cell value in patch */
template <typename T>
typename Patch<T>::CellType& Patch<T>::At(int x, int y)
{
    assert(this->IsValid());
    assert(this->IsInside(x, y));

    if (!this->mIsUnique) {
        /* If this patch is shared among several particles,
         * allocate new grid cells first because
         * this patch will be modified for the current particle and
         * will contain different cell values from other patches */
        const int numOfCells = this->mSize * this->mSize;
        CellType* pNewData = new CellType[numOfCells];
        assert(pNewData != nullptr);

        /* Copy the grid cell values to the newly allocated ones */
        std::copy_n(this->mData.get(), numOfCells, pNewData);

        /* Replace the grid cells with the newly created ones */
        this->mData.reset(pNewData);

        /* The grid cells in this patch is now unique
         * and is not shared by any other particles */
        this->mIsUnique = true;
    }

    return this->mData[y * this->mSize + x];
}

/* Get the cell value in patch */
template <typename T>
const typename Patch<T>::CellType& Patch<T>::At(int x, int y) const
{
    assert(this->IsValid());
    assert(this->IsInside(x, y));

    return this->mData[y * this->mSize + x];
}

/* Get the cell value of the patch */
template <typename T>
const typename Patch<T>::ValueType Patch<T>::Value(int x, int y) const
{
    assert(this->IsValid());
    assert(this->IsInside(x, y));
    
    return this->mData[y * this->mSize + x].Value();
}

/* Get the cell value of the patch
 * Return the default value if not allocated */
template <typename T>
const typename Patch<T>::ValueType Patch<T>::Value(
    int x, int y, const ValueType defaultVal) const
{
    /* Return the default occupancy value if not allocated */
    if (!this->IsValid())
        return defaultVal;

    assert(this->IsInside(x, y));
    return this->mData[y * this->mSize + x].Value();
}

/* Get the internal value of the specified grid cell */
template <typename T>
typename Patch<T>::StorageType Patch<T>::RawValue(
    const int x, const int y) const
{
    assert(this->IsValid());
    assert(this->IsInside(x, y));

    return this->mData[y * this->mSize + x].RawValue();
}

/* Get the internal value of the specified grid cell */
template <typename T>
typename Patch<T>::StorageType Patch<T>::RawValue(
    const int x, const int y, const StorageType defaultVal) const
{
    /* Return the default value if not allocated */
    if (!this->IsValid())
        return defaultVal;

    assert(this->IsInside(x, y));
    return this->mData[y * this->mSize + x].RawValue();
}

/* Reset all cells inside the patch */
template <typename T>
void Patch<T>::Reset()
{
    /* Do not reset if not allocated */
    if (!this->IsValid())
        return;

    /* Reset all cells inside the patch */
    for (int idxY = 0; idxY < this->mSize; ++idxY)
        for (int idxX = 0; idxX < this->mSize; ++idxX)
            this->At(idxX, idxY).Reset();
}

/* Validate the patch */
template <typename T>
bool Patch<T>::IsValid() const
{
    return (this->mSize > 0) && (this->mData != nullptr);
}

/* Check if the index is inside the patch */
template <typename T>
bool Patch<T>::IsInside(int cellX, int cellY) const
{
    return (cellX >= 0 && cellX < this->mSize) &&
           (cellY >= 0 && cellY < this->mSize);
}

} /* namespace MyGMapping */

#endif /* MY_GMAPPING_GRID_MAP_GRID_MAP_PATCH_HPP */
