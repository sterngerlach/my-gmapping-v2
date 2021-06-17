
/* gnuplot_helper.cpp */

#include "my_gmapping/io/gnuplot_helper.hpp"

#include <algorithm>
#include <cassert>

namespace MyGMapping {
namespace IO {

/* Default constructor */
GnuplotHelper::GnuplotHelper() :
    mGnuplot(nullptr)
{
    /* Interprocess communication between Gnuplot */
    this->mGnuplot.reset(popen("gnuplot 2>/dev/null", "w"));
    assert(this->mGnuplot != nullptr);

    std::fprintf(this->mGnuplot.get(), "set terminal qt\n");
    std::fprintf(this->mGnuplot.get(), "unset key\n");
    std::fprintf(this->mGnuplot.get(), "set size ratio -1\n");
}

/* Draw the particle trajectory */
void GnuplotHelper::DrawParticleTrajectory(
    const std::vector<Mapping::TimeStampedPose>& trajectory) const
{
    /* Setup trajectory nodes and edges */
    std::fprintf(this->mGnuplot.get(), "$trajectoryEdges << EOF\n");

    auto prevPoseIt = trajectory.cbegin();
    auto poseIt = std::next(prevPoseIt, 1);

    for (; poseIt != trajectory.cend(); ++prevPoseIt, ++poseIt)
        std::fprintf(this->mGnuplot.get(), "%f %f\n%f %f\n\n",
                     prevPoseIt->mPose.mX, prevPoseIt->mPose.mY,
                     poseIt->mPose.mX, poseIt->mPose.mY);

    std::fprintf(this->mGnuplot.get(), "EOF\n");

    /* Draw the trajectory of the best particle */
    std::fprintf(
        this->mGnuplot.get(),
        "plot $trajectoryEdges using 1:2 "
        "with lines lc rgb \"black\" lw 2, \\\n"
        "$trajectoryEdges using 1:2:(0.15) "
        "with circles fill solid lc rgb \"red\"\n");
    std::fflush(this->mGnuplot.get());
}

} /* namespace IO */
} /* namespace MyGMapping */
