
/* gnuplot_helper.hpp */

#ifndef MY_GMAPPING_IO_GNUPLOT_HELPER_HPP
#define MY_GMAPPING_IO_GNUPLOT_HELPER_HPP

#include <cstdio>
#include <memory>
#include <vector>

#include "my_gmapping/pose.hpp"
#include "my_gmapping/mapping/trajectory_node.hpp"

namespace MyGMapping {
namespace IO {

/*
 * PipeDeleter struct is a custom deleter for std::unique_ptr
 */
struct PipeDeleter
{
    void operator()(FILE* ptr) const
    {
        if (ptr != nullptr)
            pclose(ptr);
    }
};

/*
 * GnuplotHelper class is for visualizing pose graphs to the Gnuplot window
 */
class GnuplotHelper
{
public:
    /* Default constructor */
    GnuplotHelper();
    /* Destructor */
    ~GnuplotHelper() = default;

    /* Copy constructor (disabled) */
    GnuplotHelper(const GnuplotHelper&) = delete;
    /* Copy assignment operator (disabled) */
    GnuplotHelper& operator=(const GnuplotHelper&) = delete;
    /* Move constructor (disabled) */
    GnuplotHelper(GnuplotHelper&&) = delete;
    /* Move assignment operator (disabled) */
    GnuplotHelper& operator=(GnuplotHelper&&) = delete;

    /* Draw the particle trajectory */
    void DrawParticleTrajectory(
        const std::vector<Mapping::TimeStampedPose>& trajectory) const;

private:
    std::unique_ptr<FILE, PipeDeleter> mGnuplot;
};

} /* namespace IO */
} /* namespace MyGMapping */

#endif /* MY_GMAPPING_IO_GNUPLOT_HELPER_HPP */
