
/* bresenham.hpp */

#ifndef MY_GMAPPING_BRESENHAM_HPP
#define MY_GMAPPING_BRESENHAM_HPP

#include <cmath>
#include <vector>

#include "my_gmapping/point.hpp"

namespace MyGMapping {

/* Perform the Bresenham algorithm */
void Bresenham(const Point2D<int>& startIdx,
               const Point2D<int>& endIdx,
               std::vector<Point2D<int>>& indices);

/* Perform the Bresenham algorithm at subpixel accuracy */
void BresenhamScaled(const Point2D<int>& scaledStartIdx,
                     const Point2D<int>& scaledEndIdx,
                     const int subpixelScale,
                     std::vector<Point2D<int>>& indices);

} /* namespace MyGMapping */

#endif /* MY_GMAPPING_BRESENHAM_HPP */
