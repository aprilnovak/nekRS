/**
 * Structure containing boundary condition information that is used to
 * apply boundary conditions in your self-defined functions in the
 * .oudf file.
 */
struct bcData
{
  int idM;
  int fieldOffset;
  int id;

  int scalarId;

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

  dfloat uM, vM, wM;
  dfloat uP, vP, wP;
  dfloat uxP, uyP, uzP;
  dfloat vxP, vyP, vzP;
  dfloat wxP, wyP, wzP;

  dfloat pM;
  dfloat pP;

  @globalPtr const dfloat* wrk;

  dfloat sM, sP, sF;
};
