/*
 *
 * A model to impose weak periodic boundary 
 * conditions on the 2D unit cells.
 *
 * Erik Giesen Loo, April 2018
 *
 */

#ifndef WEAKPERIODICBC_MODEL_H
#define WEAKPERIODICBC_MODEL_H

//-----------------------------------------------------------------------
//   class WeakPeriodicBCModel
//-----------------------------------------------------------------------

#include <jive/model/Model.h>
#include <jive/model/ModelFactory.h>

#include <jive/fem/NodeGroup.h>
#include <jive/fem/XNodeSet.h>
#include <jive/util/XDofSpace.h>
#include <jive/util/Constraints.h>
#include <jive/util/Assignable.h>

#include <jive/algebra/MatrixBuilder.h>
#include <jive/geom/BoundaryLine.h>
#include <jive/geom/Line.h>

using namespace jem;

using jem::numeric::Function;
using jem::util::Properties;
using jive::BoolVector;
using jive::IdxVector;
using jive::Vector;
using jive::algebra::MatrixBuilder;
using jive::fem::NodeSet;
using jive::fem::XNodeSet;
using jive::geom::BoundaryLine2;
using jive::geom::BoundaryShape;
using jive::geom::Line;
using jive::model::Model;
using jive::util::Assignable;
using jive::util::Constraints;
using jive::util::XDofSpace;

//=========================================================
//    class WeakPeriodicBCModel
//=========================================================

class WeakPeriodicBCModel : public Model
{
public: // public members:
  typedef WeakPeriodicBCModel Self;
  typedef Model Super;
  typedef Array<Ref<Function>, 1> FuncVector;

  static const char *NODE_GROUPS[6];
  static const char *STRAINRATE_PROP;
  static const char *STRAINPATH_PROP;
  static const char *MAXTIME_PROP;
  static const char *ACTIVE_PROP;

  enum StrainType
  {
    Free,
    ConstantRate,
    Path
  };

public: // public methods:
  WeakPeriodicBCModel(const String &name,
                      const Properties &conf,
                      const Properties &props,
                      const Properties &globdat);

  virtual void configure(const Properties &props,
                         const Properties &globdat);

  virtual void getConfig(const Properties &conf,
                         const Properties &globdat) const;

  virtual bool takeAction(const String &action,
                          const Properties &params,
                          const Properties &globdat);

protected: // Protected Methods (internal use)
  virtual ~WeakPeriodicBCModel();

private: // Private Methods (internal use)
  void init_(const Properties &globdat);
  void advance_() const;
  void checkCommit_(const Properties &params, const Properties &globdat) const;
  void fixCorner_() const;
  void applyStrain_(const Vector &strain) const;
  void setPeriodicCons_() const;
  FuncVector makeStrainFuncs_(const Vector &strainRate) const;
  FuncVector getStrainFuncs_(const Properties &globdat) const;
  void sortBndNodes_();
  void sortBndFace_(IdxVector &bndFace, const idx_t &index);
  void createTractionMesh_();
  void coarsenMesh_(IdxVector &trFace, const idx_t &index);
  void augmentMatrix_(Ref<MatrixBuilder> mbuilder,
                      const Vector &force,
                      const Vector &disp);
  void getTractionMeshNodes_(IdxVector &connect,
                             const Vector &x,
                             const idx_t &face);

private: // Private Members (internal use)
  Assignable<XNodeSet> nodes_;
  BoolVector active_;
  idx_t rank_; // number of dimensions

  Ref<XDofSpace> dofs_;
  Ref<Constraints> cons_;
  IdxVector U_doftypes_;      // displcement dof types
  IdxVector T_doftypes_;      // traction dof types
  Ref<BoundaryShape> bshape_; // boundary element
  int nIP_;                   // number of int. points of boundary element
  int nnod_;                  // number of nodes of boundary element
  int localrank_;             // local rank of boundary element

  IdxVector bndNodes_[6];   // boundary nodes [ xmin, xmax, ymin, ymax ]
  IdxVector trNodes_[3];    // traction mesh nodes [ xmin, ymin ]
  Tuple<idx_t, 3> masters_; // corner nodes [ cornerX, cornerY, cornerZ ]
  idx_t ifixed_;            // master corner node [ corner0 ]

  Vector imposedStrain_; // total applied strain

  double time_;
  double stepSize_;
  double maxTime_;

  Vector box_;   // specimen coordinates
  Vector dx_;    // specimen dimensions
  Vector dx0_;   // smallest element dimensions
  double factor; // coarsening factor

  String strainFile_;
  StrainType strainType_;
  FuncVector strainFunc_;
};

#endif