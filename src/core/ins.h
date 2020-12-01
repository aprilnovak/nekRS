#if !defined(nekrs_ins_hpp_)
#define nekrs_ins_hpp_

#include "mesh2D.h"
#include "mesh3D.h"
#include "elliptic.h"
#include "cds.h"
#include "linAlg.hpp"

/**
 * Structure that holds most of the high-level information regarding
 * the equation set to be solved, the mesh, and the solution method.
 */
typedef struct
{
  /// Mesh dimension
  int dim;

  /// Element type
  int elementType;

  mesh_t* mesh;
  mesh_t* meshT;

  elliptic_t* uSolver;
  elliptic_t* vSolver;
  elliptic_t* wSolver;
  elliptic_t* uvwSolver;
  elliptic_t* pSolver;

  /**
   * \brief Convection-diffusion solver for passive scalars
   *
   * If no passive scalars are present in the simulation, this will
   * not be initialized.
   */
  cds_t* cds;

  oogs_t* gsh;

  linAlg_t* linAlg;

  dlong ellipticWrkOffset;

  /// Whether the equation set includes flow, i.e. a solution for velocity
  int flow;

  /// Number of scalars to solve for, or the number of passive scalars plus temperature
  int Nscalar;

  setupAide options;

  setupAide vOptions;

  /// User-specified options for the pressure solver
  setupAide pOptions;

  /// Number of velocity fields solved for, which is equal to the mesh dimension, \ref ins_t.dim
  int NVfields;

  /**
   * \brief Total number of flow-related fields to solve for
   *
   * This is computed as the number of velocity components, \ref ins_t.NVfields, plus pressure
   */
  int NTfields;

  /**
   * \brief Size of a solution field
   *
   * Size of a solution field (a velocity component, pressure, or a passive scalar) in terms of
   * the number of floats required to represent it. This is taken as the maximum size needed
   * to represent any of the various solution fields in the case that conjugate heat transfer
   * exists (and you therefore have different meshes for the fluid and solid phases). This
   * value is adjusted so that it is an even number of "pages", where the page size is
   * by default assumed to be 4096, but can also be specified with the NEKRS_PAGE_SIZE
   * environment variable.
   */
  dlong fieldOffset;

  /// Number of GLL points local to this process
  dlong Nlocal;

  dlong Ntotal;

  int Nblock;


  /// Time step size
  dfloat dt;

  /// Storage of \f$\frac{1}{\Delta t}\f$ for fast evaluation in kernels,
  /// where \f$\Delta t\f$ is the time step
  dfloat idt;
  dfloat time;
  int tstep;
  dfloat g0, ig0;

  /// Time at which to start the simulation
  dfloat startTime;

  /** 
   * \brief Time at which to end the simulation
   *
   * If \ref ins_t.startTime is not zero, this
   * is adjusted to the sum of the user-specified final time and the start time.
   * Note also that if the difference between the final time and the start time
   * is not evenly divisible by the requested time step, that one extra time
   * step is added such that the simulation end time will be slightly larger than
   * the specified final time.
   */
  dfloat finalTime;

  /// Whether the equation set includes conjugate heat transfer
  int cht;

  /// Order of accuracy of the time integrator
  int temporalOrder;

  int ExplicitOrder;

  /**
   * \brief Total number of time steps in the simulation
   *
   * The total number of time steps is evaluated based on the
   * start time, final time, and time step size. If the difference between the
   * final time and start time is not evenly divisible by the time step size,
   * one extra time step is added to ensure simulation up to and including the
   * final time.
   */
  int NtimeSteps;

  /// Number of stages for the time integrator
  int Nstages;

  /// Time step interval on which to output the simulation
  int outputStep;

  /// Whether simulation data is to be written on this time step
  int isOutputStep;

  int outputForceStep;

  /// Number of iterations performed for most recent \f$x\f$-velocity solve
  int NiterU;

  /// Number of iterations performed for most recent \f$y\f$-velocity solve
  int NiterV;

  /// Number of iterations performed for most recent \f$z\f$-velocity solve
  int NiterW;

  /// Number of iterations performed for most recent pressure solve
  int NiterP;
  dfloat presTOL, velTOL;

  ///@{ 
  /**
   * \brief Velocity solution for all components on the host and device
   *
   * This holds the velocity solution for all components in the \ref ins_t.dim
   * space, for each of the \ref ins_t.Nstages of the time integrator (that is,
   * this array holds all the previous time step information needed for the
   * time integration (but not necessarily _all_ previous time steps - only what
   * is needed to move forward in time).
   */
  dfloat* U;
  occa::memory o_U;
  ///@}

  ///@{
  /**
   * \brief Pressure solution on the host and device
   *
   * This holds the pressure solution for each of the \ref ins_t.Nstages of the time
   * integrator (that is, this array holds all the previous time step information
   * needed for the time integration (but not necessarily _all_ previous time steps -
   * only what is needed to move forward in time).
   */
  dfloat * P;
  occa::memory o_P;
  ///@}

  dfloat* BF, * FU;

  //RK Subcycle Data
  int SNrk;
  dfloat* Srka, * Srkb, * Srkc;
  occa::memory o_Srka, o_Srkb;

  //ARK data
  int Nrk;
  dfloat* rkC;

  //EXTBDF data
  dfloat* extbdfA, * extbdfB, * extbdfC;
  dfloat* extC;

  int* VmapB;
  occa::memory o_VmapB;

  dlong* elementInfo;
  occa::memory o_elementInfo;

  occa::memory o_wrk0, o_wrk1, o_wrk2, o_wrk3, o_wrk4, o_wrk5, o_wrk6, o_wrk7,
               o_wrk9, o_wrk12, o_wrk15;

  /// Number of subcycling time steps
  int Nsubsteps;
  dfloat* Ue, sdt;
  occa::memory o_Ue;

  dfloat* div;
  occa::memory o_div;

  /// Fluid density
  dfloat rho;

  /// Fluid viscosity
  dfloat mue;
  occa::memory o_rho, o_mue;

  //@{
  /**
   * \brief Arrays in user space on the host and device
   *
   * These arrays are free for the use of the user - no other routines touch
   * these arrays. Because these arrays are free for whatever purpose the user
   * desires, the convention is to allocate space for these arrays and write to
   * them in the user-defined function (.udf) file. To have access to these arrays
   * for the entire duration of the simulation, you can allocate space for them
   * in `UDF_Setup` and write to them for each time step in `UDF_ExecuteStep`.
   */
  dfloat* usrwrk;
  occa::memory o_usrwrk;
  //@}

  occa::memory o_idH; // i.e. inverse of 1D Gll Spacing for quad and Hex

  /// Whether to read from a restart file
  int readRestartFile;

  /// Whether to write a restart file
  int writeRestartFile;

  int restartedFromFile;

  int filterNc; // filter cut modes i.e. below is not touched
  dfloat* filterM, filterS;
  occa::memory o_filterMT; // transpose of filter matrix
  occa::kernel filterRTKernel; // Relaxation-Term based filtering

  occa::kernel qtlKernel;
  occa::kernel pressureAddQtlKernel;
  occa::kernel pressureStressKernel;

  occa::kernel PQKernel;
  occa::kernel mueDivKernel;
  occa::kernel dotMultiplyKernel;

  occa::kernel scalarScaledAddKernel;
  occa::kernel scaledAddKernel;
  occa::kernel subCycleVolumeKernel,  subCycleCubatureVolumeKernel;
  occa::kernel subCycleSurfaceKernel, subCycleCubatureSurfaceKernel;
  occa::kernel subCycleRKUpdateKernel;
  occa::kernel extrapolateKernel;

  occa::kernel wgradientVolumeKernel;

  occa::kernel subCycleStrongCubatureVolumeKernel;
  occa::kernel subCycleStrongVolumeKernel;

  occa::kernel constrainKernel;

  occa::memory o_BF;
  occa::memory o_FU;

  int var_coeff;
  dfloat* prop, * ellipticCoeff;
  occa::memory o_prop, o_ellipticCoeff;

  occa::memory o_UH;

  occa::memory o_vHaloBuffer, o_pHaloBuffer;
  occa::memory o_velocityHaloGatherTmp;

  occa::kernel haloGetKernel;
  occa::kernel haloPutKernel;

  //ARK data
  occa::memory o_rkC;

  //EXTBDF data
  occa::memory o_extbdfA, o_extbdfB, o_extbdfC;
  occa::memory o_extC;

  occa::kernel advectionVolumeKernel;
  occa::kernel advectionCubatureVolumeKernel;

  occa::kernel advectionStrongVolumeKernel;
  occa::kernel advectionStrongCubatureVolumeKernel;

  occa::kernel diffusionKernel;
  occa::kernel velocityGradientKernel;

  occa::kernel gradientVolumeKernel;

  occa::kernel divergenceVolumeKernel;
  occa::kernel divergenceSurfaceKernel;

  occa::kernel divergenceStrongVolumeKernel;
  occa::kernel sumMakefKernel;
  occa::kernel pressureRhsKernel;
  occa::kernel pressureDirichletBCKernel;
  occa::kernel pressurePenaltyKernel;
  occa::kernel pressureUpdateKernel;

  occa::kernel velocityRhsKernel;
  occa::kernel velocityNeumannBCKernel;
  occa::kernel velocityDirichletBCKernel;

  occa::kernel fillKernel;

  occa::kernel cflKernel;
  occa::kernel maxKernel;

  occa::kernel setEllipticCoeffKernel;
  occa::kernel setEllipticCoeffPressureKernel;

  occa::kernel pressureAxKernel;
  occa::kernel curlKernel;
  occa::kernel invMassMatrixKernel;
  occa::kernel massMatrixKernel;

  occa::kernel maskCopyKernel;

  int* EToB;
  occa::memory o_EToB;

  occa::properties* kernelInfo;
}ins_t;

#endif
