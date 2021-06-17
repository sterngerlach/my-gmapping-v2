
/* grid_map_base.hpp */

#ifndef MY_GMAPPING_GRID_MAP_GRID_MAP_BASE_HPP
#define MY_GMAPPING_GRID_MAP_GRID_MAP_BASE_HPP

#include "my_gmapping/point.hpp"
#include "my_gmapping/pose.hpp"

namespace MyGMapping {

/*
 * GridMapBase is the base class for the occupancy grid maps and provides
 * useful methods for reading the grid cells, accessing the map parameters,
 * and transforming coordinates, however it does not have any methods for
 * updating the grid cells nor accessing the patches (small chunk of the map),
 * thus the grid map is read-only
 */
template <typename T, typename U>
class GridMapBase
{
public:
    /* Constructor */
    GridMapBase() = default;
    /* Destructor */
    virtual ~GridMapBase() = default;

    /* Copy constructor */
    GridMapBase(const GridMapBase&) = default;
    /* Copy assignment operator */
    GridMapBase& operator=(const GridMapBase&) = default;
    /* Move constructor */
    GridMapBase(GridMapBase&&) noexcept = default;
    /* Move assignment operator */
    GridMapBase& operator=(GridMapBase&&) noexcept = default;

    /* Get the unknown occupancy probability value */
    virtual T UnknownValue() const = 0;
    /* Get the internal value to represent unknown occupancy probability */
    virtual U UnknownRawValue() const = 0;
    /* Get the minimum internal value */
    virtual U RawValueMin() const = 0;
    /* Get the maximum internal value */
    virtual U RawValueMax() const = 0;
    /* Get the internal value range */
    virtual U RawValueRange() const = 0;

    /* Check if the index is inside the map */
    virtual bool IsInside(int cellX, int cellY) const = 0;
    /* Check if the index is inside the map */
    virtual bool IsInside(const Point2D<int>& cellIdx) const = 0;

    /* Check if the point is inside the map */
    virtual bool IsInside(double mapX, double mapY) const = 0;
    /* Check if the point is inside the map */
    virtual bool IsInside(const Point2D<double>& mapPos) const = 0;

    /* Check if the cell is allocated */
    virtual bool IsAllocated(int cellX, int cellY) const = 0;
    /* Check if the cell is allocated */
    virtual bool IsAllocated(const Point2D<int>& cellIdx) const = 0;

    /* Convert cell index to map coordinate
     * Calculate the minimum position of the specified cell (x, y) */
    virtual Point2D<double> CellIndexToMapCoordinate(int x, int y) const = 0;
    /* Convert cell index to map coordinate
     * Calculate the minimum position of the specified cell (x, y) */
    virtual Point2D<double> CellIndexToMapCoordinate(
        const Point2D<int>& cellIdx) const = 0;

    /* Convert map coordinate to cell index */
    virtual Point2D<int> MapCoordinateToCellIndex(
        const double x, const double y) const = 0;
    /* Convert map coordinate to cell index */
    virtual Point2D<int> MapCoordinateToCellIndex(
        const Point2D<double>& mapPos) const = 0;

    /* Convert map coordinate to floating-point cell index */
    virtual Point2D<double> MapCoordinateToCellIndexFloat(
        const double x, const double y) const = 0;
    /* Convert map coordinate to floating-point cell index */
    virtual Point2D<double> MapCoordinateToCellIndexFloat(
        const Point2D<double>& mapPos) const = 0;

    /* Get the value of the specified index */
    virtual T Value(int x, int y) const = 0;
    virtual T Value(const Point2D<int>& cellIdx) const = 0;

    /* Get the value of the specified index */
    virtual T Value(int x, int y, const T defaultVal) const = 0;
    virtual T Value(const Point2D<int>& cellIdx,
                    const T defaultVal) const = 0;

    /* Get the internal value of the specified grid cell */
    virtual U RawValue(const int idxX, const int idxY) const = 0;
    /* Get the internal value of the specified grid cell */
    virtual U RawValue(const Point2D<int>& cellIdx) const = 0;

    /* Get the internal value of the specified grid cell
     * The default value is returned if the specified grid cell is
     * out of bounds or is not yet allocated */
    virtual U RawValue(const int idxX, const int idxY,
                       const U defaultVal) const = 0;
    /* Get the internal value of the specified grid cell */
    virtual U RawValue(const Point2D<int>& cellIdx,
                       const U defaultVal) const = 0;

    /* Calculate the distance between two cell indices */
    virtual double Distance(int x0, int y0, int x1, int y1) const = 0;
    virtual double Distance(const Point2D<int>& cellIdx0,
                            const Point2D<int>& cellIdx1) const = 0;

    /* Calculate the squared distance between two cell indices */
    virtual double SquaredDistance(int x0, int y0, int x1, int y1) const = 0;
    virtual double SquaredDistance(const Point2D<int>& cellIdx0,
                                   const Point2D<int>& cellIdx1) const = 0;

    /* Get the cell size */
    virtual double CellSize() const = 0;

    /* Get the number of the cells */
    virtual int NumCellsX() const = 0;
    virtual int NumCellsY() const = 0;

    /* Get the size of the map */
    virtual double MapSizeX() const = 0;
    virtual double MapSizeY() const = 0;

    /* Get the minimum position of the grid map (in world coordinate) */
    virtual const Point2D<double>& MinPos() const = 0;
};

} /* namespace MyGMapping */

#endif /* MY_GMAPPING_GRID_MAP_GRID_MAP_BASE_HPP */
