#if !defined(nekrs_udf_hpp_)
#define nekrs_udf_hpp_

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

extern "C" {

/**
 * \brief User-defined function called immediately after the JIT compilation of the user-defined functions.
 *
 * This calls the \ref UDF.setup0 function pointer. Because this is a user-defined
 * function, you can use this function for whatever you like, but it is typically
 * used to set options in the `setupAide` for continuous integration testing where
 * the problem options may be changed from what is set in the `.par` file. To define
 * this function, add it to your case's `.udf` file.
 *
 * @param[in] comm         MPI communicator
 * @param[in, out] options options to set up the simulation
 **/
void UDF_Setup0(MPI_Comm comm, setupAide &options);

/**
 * \brief User-defined function called once at the beginning of the simulation.
 *
 * This calls the \ref UDF.setup function pointer. And more precisely, this function
 * is called just after initializing the flow solver and reading the restart file (if
 * present). Because this is a user-defined function, you can use this function for
 * whatever you like, but it is typically used to set up scratch space (the
 * `usrwrk` arrays) and get an initial condition from Nek5000's `useric` routine (which
 * would be in a Fortran `.usr` file) if a restart file is not available. You could also
 * write initial conditions directly in this function if you prefer to stay with C++.
 *
 * In this routine, we can also set the other function pointers in \ref UDF that aren't
 * set in `udfLoad` (i.e. \ref UDF.uEqnSource, \ref UDF.sEqnSource, \ref UDF.properties,
 * and \ref UDF.div.
 * @param[in] ins flow solver
 **/
void UDF_Setup(ins_t* ins);

void UDF_LoadKernels(ins_t* ins);

/**
 * \brief User-defined function called once at the beginning of the simulation and then once per step
 *
 * When called at the beginning of the simulation, this is called just before copying
 * the solution on the device to the host, and is the last action that happens in the
 * simulation setup phase. Within each time step, this is called after solving for
 * that time step.
 * @param[in] ins   flow solver
 * @param[in] time  current simulation time
 * @param[in] tstep time step number
 **/
void UDF_ExecuteStep(ins_t* ins, dfloat time, int tstep);
};

/**
 * Function pointer to a function taking two inputs and returning void
 * @param[in] comm         MPI communicator
 * @param[in, out] options options to set up the simulation
 **/
typedef void (* udfsetup0)(MPI_Comm comm, setupAide &options);

/**
 * Function pointer to a function taking one input and returning void
 * @param[in] ins flow solver
 **/
typedef void (* udfsetup)(ins_t* ins);

/**
 * Function pointer to a function taking one input and returning void
 * @param[in] ins flow solver
 **/
typedef void (* udfloadKernels)(ins_t* ins);

/**
 * Function pointer to a function taking three inputs and returning void
 * @param[in] ins   flow solver
 * @param[in] time  current solution time
 * @param[in] tstep time step number
 **/
typedef void (* udfexecuteStep)(ins_t* ins, dfloat time, int tstep);

typedef void (* udfuEqnSource)(ins_t* ins, dfloat time, occa::memory o_U, occa::memory o_FU);
typedef void (* udfsEqnSource)(ins_t* ins, dfloat time, occa::memory o_S, occa::memory o_SU);
typedef void (* udfproperties)(ins_t* ins, dfloat time, occa::memory o_U,
                               occa::memory o_S, occa::memory o_UProp,
                               occa::memory o_SProp);
typedef void (* udfdiv)(ins_t* ins, dfloat time, occa::memory o_div);

/// Function pointers to a number of user-defined functions
typedef struct
{
  /// Function pointer to `UDF_Setup0`
  udfsetup0 setup0;

  /// Function pointer to `UDF_Setup`
  udfsetup setup;

  /// Function pointer to `UDF_LoadKernels`
  udfloadKernels loadKernels;

  /// Function pointer to `UDF_ExecuteStep`
  udfexecuteStep executeStep;

  /**
   * \brief Function pointer to any function that gives the momentum equation source
   *
   * Note that the function must have the correct input/output parameter signature to
   * match the type definition for \ref udfuEqnSource.
   **/
  udfuEqnSource uEqnSource;

  /**
   * \brief Function pointer to any function that gives the scalar equation source
   *
   * Note that the function must have the correct input/output parameter signature to
   * match the type definition for \ref udfsEqnSource.
   **/
  udfsEqnSource sEqnSource;

  /**
   * \brief Function pointer to any function that gives the material properties
   *
   * Note that the function must have the correct input/output parameter signature to
   * match the type definition for \ref udfproperties.
   **/
  udfproperties properties;

  udfdiv div;
} UDF;

extern UDF udf;

void udfBuild(const char* udfFile);

/**
 * \brief Assign the function pointers to four of the user-defined functions
 *
 * This function assigns the \ref UDF.setup0, \ref UDF.setup, \ref UDF.loadKernels,
 * and \ref UDF.executeStep function pointers. Based on the error checking flags
 * passed to the `udfLoadFunction` routine, only the `UDF_Setup0` function is optional.
 * Note that this function is only called if a `.udf` file is present for your case.
 **/
void udfLoad(void);

/**
 * Return a function pointer to a specified function
 * @param[in] fname  Name of the function
 * @param[in] errchk Whether to print an error if function with `fname` can't be found
 **/
void* udfLoadFunction(const char* fname, int errchk);

occa::kernel udfBuildKernel(ins_t* ins, const char* function);

#endif
