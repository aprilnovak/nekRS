#ifndef CDS_H
#define CDS_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"
#include "mesh3D.h"
#include "elliptic.h"

#define NSCALAR_MAX 100

/**
 * Structure containing most of the information related to the convection-diffusion
 * passive scalar solve. Most of these fields are directly copied from those in the
 * \ref ins_t solver.
 */
typedef struct
{
  /// Mesh dimension
  int dim, elementType;

  mesh_t* mesh;

  mesh_t* meshV;

  elliptic_t* solver[NSCALAR_MAX];

  /// Number of velocity fields
  int NVfields;

  /// Number of scalars to solve for, or the number of passive scalars plus temperature
  int NSfields;

  setupAide options;

  oogs_t *gsh, *gshT;

  dlong vFieldOffset;
  dlong fieldOffset;
  dlong Nlocal, Ntotal;
  int Nblock;
  dfloat dt, idt, cfl, dti;          // time step
  dfloat time;
  int tstep;
  dfloat g0, ig0;
  dfloat startTime;
  dfloat finalTime;

  int temporalOrder;
  int ExplicitOrder;
  int NtimeSteps;  // number of time steps
  int Nstages;
  int outputStep;
  int outputForceStep;
  int dtAdaptStep;

  int compute[NSCALAR_MAX];

  /// Number of iterations taken in most recent scalar solve for each passive scalar
  int Niter[NSCALAR_MAX];

  //solver tolerances
  dfloat TOL;

  /**
   * \brief Velocity solution for all components
   *
   * This is simply a pointer to \ref ins_t.U, which is used to transport the passive scalars.
   */
  dfloat* U;

  ///@{
  /**
   * \brief Scalar solution fields on the host and device
   *
   * These hold the solution for all \ref cds_t.NSfields passive scalars, for each
   * of the \ref cds_t.Nstages of the time integrator (that is, these arrays hold all
   * the previous time step information needed for the time integration (but not
   * necessarily _all_ previous time steps - only what is needed to move forward in time.)
   */
  dfloat* S;
  occa::memory o_S;
  ///@}

  dfloat* rkNS;
  //  dfloat *rhsS;
  dfloat* rkS;

  //RK Subcycle Data
  int SNrk;
  dfloat* Srka, * Srkb, * Srkc;
  occa::memory o_Srka, o_Srkb;

  //EXTBDF data
  dfloat* extbdfA, * extbdfB, * extbdfC;
  dfloat* extC;

  int* mapB[NSCALAR_MAX], * EToB[NSCALAR_MAX];
  occa::memory o_mapB[NSCALAR_MAX];
  occa::memory o_EToB[NSCALAR_MAX];

  /**
   * \brief Array in user space on the device
   *
   * This array is free for the use of the user - no other routines touch
   * this array.
   */
  occa::memory* o_usrwrk;

  //halo data
  dfloat* sendBuffer;
  dfloat* recvBuffer;
  dfloat* haloGatherTmp;
  // //
  dfloat* ssendBuffer;
  dfloat* srecvBuffer;
  dfloat* shaloGatherTmp;

  occa::memory o_sendBuffer, h_sendBuffer;
  occa::memory o_recvBuffer, h_recvBuffer;
  occa::memory o_gatherTmpPinned, h_gatherTmpPinned;

  //
  occa::memory o_ssendBuffer, h_ssendBuffer;
  occa::memory o_srecvBuffer, h_srecvBuffer;
  occa::memory o_sgatherTmpPinned, h_sgatherTmpPinned;

  int Nsubsteps;
  dfloat sdt;
  dfloat* Ue;
  occa::memory o_Ue;

  int var_coeff;
  dfloat* prop, * ellipticCoeff;
  occa::memory o_prop, o_ellipticCoeff;
  occa::memory o_rho, o_diff;

  dfloat* cU, * cSd, * cS, * FS, * BF;
  occa::memory o_cU, o_cSd, o_cS, o_FS, o_BF, o_BFDiag;

  occa::memory o_wrk0, o_wrk1, o_wrk2, o_wrk3, o_wrk4, o_wrk5, o_wrk6;

  occa::kernel fillKernel;
  occa::kernel sumMakefKernel;
  occa::kernel scaledAddKernel;
  occa::kernel subCycleVolumeKernel,  subCycleCubatureVolumeKernel;
  occa::kernel subCycleSurfaceKernel, subCycleCubatureSurfaceKernel;
  occa::kernel subCycleRKUpdateKernel;
  occa::kernel subCycleStrongCubatureVolumeKernel;
  occa::kernel subCycleStrongVolumeKernel;

  occa::kernel filterRTKernel; // Relaxation-Term based filtering
  // occa::kernel constrainKernel;

  occa::memory o_U;
  occa::memory o_Se;

  // occa::memory o_Vort, o_Div; // Not sure to keep it
  occa::memory o_haloBuffer;
  occa::memory o_haloGatherTmp;

  occa::memory o_shaloBuffer;
  occa::memory o_shaloGatherTmp;

  //ARK data
  occa::memory o_rkC;

  //EXTBDF data
  occa::memory o_extbdfA, o_extbdfB, o_extbdfC;
  occa::memory o_extC;

// Will be depreceated.....AK
  occa::kernel haloExtractKernel;
  occa::kernel haloScatterKernel;
  occa::kernel scalarHaloExtractKernel;
  occa::kernel scalarHaloScatterKernel;

  occa::kernel haloGetKernel;
  occa::kernel haloPutKernel;
  occa::kernel scalarHaloGetKernel;
  occa::kernel scalarHaloPutKernel;

  occa::kernel setFlowFieldKernel;
  occa::kernel setScalarFieldKernel;

  occa::kernel advectionVolumeKernel;
  occa::kernel advectionSurfaceKernel;
  occa::kernel advectionCubatureVolumeKernel;
  occa::kernel advectionCubatureSurfaceKernel;
  occa::kernel advectionStrongVolumeKernel;
  occa::kernel advectionStrongCubatureVolumeKernel;

  occa::kernel helmholtzRhsIpdgBCKernel;
  occa::kernel helmholtzRhsBCKernel;
  occa::kernel dirichletBCKernel;
  occa::kernel setEllipticCoeffKernel;

  occa::kernel invMassMatrixKernel;
  occa::kernel massMatrixKernel;

  occa::kernel maskCopyKernel;

  occa::properties* kernelInfo;
}cds_t;

occa::memory cdsSolve(int i, cds_t* cds, dfloat time);

#endif
