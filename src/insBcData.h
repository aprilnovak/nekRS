/**
 * Structure containing boundary condition information that is used to
 * apply boundary conditions in your self-defined functions in the
 * .oudf file.
 */
struct bcData
{
  int idM;
  int fieldOffset;

  /// Boundary ID to which boundary condition is being applied
  int id;

  /// Scalar ID to which boundary condition is being applied
  int scalarId;

  /// Simulation time
  dfloat time;

  /// \f$x\f$ coordinate
  dfloat x;

  /// \f$y\f$ coordinate
  dfloat y;

  /// \f$z\f$ coordinate
  dfloat z;

  /// Surface normal component in \f$x\f$ direction
  dfloat nx;

  /// Surface normal component in \f$y\f$ direction
  dfloat ny;

  /// Surface normal component in \f$z\f$ direction
  dfloat nz;

  ///@{
  /**
   * \brief \f$x\f$, \f$y\f$, and \f$z\f$ values of velocity on the negative side
   *
   * This feature is present for future implementation of discontinuous Galerkin
   * capabilities, and is currently unused in the boundary conditions.
   **/
  dfloat uM;
  dfloat vM;
  dfloat wM;
  ///@}

  ///@{
  /**
   * \brief \f$x\f$, \f$y\f$, and \f$z\f$ values of velocity on the positive side
   * 
   * These values are used to obtain the values of velocity on the boundary.
   **/
  dfloat uP;
  dfloat vP;
  dfloat wP;
  ///@}

  dfloat uxP, uyP, uzP;
  dfloat vxP, vyP, vzP;
  dfloat wxP, wyP, wzP;

  /**
   * \brief Value of pressure on the negative side
   *
   * This feature is present for future implementation of discontinuous Galerkin
   * capabilities, and is currently unused in the boundary conditions.
   **/
  dfloat pM;

  /**
   * \brief Value of pressure on the positive side
   * 
   * This value is used to obtain the value of pressure on the boundary.
   **/
  dfloat pP;

  @globalPtr const dfloat* wrk;

  /**
   * \brief Value of the scalar on the negative side
   *
   * This feature is present for future implementation of discontinuous Galerkin
   * capabilities, and is currently unused in the boundary conditions.
   **/
  dfloat sM;

  /**
   * \brief Value of the scalar on the positive side
   * 
   * This value is used to obtain the value of the scalar on the boundary.
   **/
  dfloat sP;

  /// Value of the scalar flux on the boundary
  dfloat sF;
};
