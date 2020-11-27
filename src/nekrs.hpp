#if !defined(nekrs_nrs_hpp_)
#define nekrs_nrs_hpp_

#include <string>
#include <mpi.h>

/**
 * Namespace containing all of the high-level routines and information
 * needed to run nekRS; this namespace forms the nekRS API for
 * external applications.
 */
namespace nekrs
{
void setup(MPI_Comm comm, int buildOnly, int sizeTarget,
           int ciMode, std::string cacheDir, std::string setupFile,
           std::string backend, std::string deviceID);

void runStep(double time, double dt, int tstep);

/**
 * \brief Copy the velocity, pressure, and scalar solutions from device to host
 *
 * Copy the velocity, pressure, and passive scalar solutions from the device to the
 * nekRS data structures on the host. Then, because Nek5000 is still used for I/O,
 * copy the host solutions to the `nekData` structure.
 **/
void copyToNek(double time, int tstep);

void udfExecuteStep(double time, int tstep, int isOutputStep);
void nekOutfld(void);
void nekUserchk(void);
void printRuntimeStatistics(void);

const double dt(void);
const int outputStep(void);

/**
 * Get the number of time steps in the simulation
 * \return number of time steps
 */
const int NtimeSteps(void);

/**
 * Get the start time of the simulation
 * \return start time
 */
const double startTime(void);

/**
 * Get the end time of the siulation
 * \return end time
 */
const double finalTime(void);

void* nekPtr(const char* id);
}

#endif
