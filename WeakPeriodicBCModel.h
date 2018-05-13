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

#include <jem/util/Flex.h>
#include <jive/fem/NodeGroup.h>
#include <jive/fem/XNodeSet.h>
#include <jive/util/Assignable.h>
#include <jive/util/Constraints.h>
#include <jive/util/XDofSpace.h>

#include <jive/algebra/MatrixBuilder.h>
#include <jive/geom/BoundaryLine.h>

using namespace jem;

using jem::numeric::Function;
using jem::util::Properties;
using jive::model::Model;

using jive::Vector;
using jive::IdxVector;
using jive::BoolVector;
using jem::util::Flex;
using jive::fem::NodeSet;
using jive::fem::XNodeSet;
using jive::util::Assignable;
using jive::util::Constraints;
using jive::util::XDofSpace;

using jive::algebra::MatrixBuilder;
using jive::geom::BoundaryLine2;
using jive::geom::BoundaryShape;

//=========================================================
//    class WeakPeriodicBCModel
//=========================================================

class WeakPeriodicBCModel : public Model
{
public: // Public Members:
  typedef WeakPeriodicBCModel Self;
  typedef Model Super;
  typedef Flex<idx_t> FlexVector;
  typedef Flex<idx_t>::Iterator Iter;
  typedef Array<Ref<Function>, 1> FuncVector;

  static const char *NODE_GROUPS[6];
  static const char *STRAINRATE_PROP;
  static const char *STRAINPATH_PROP;
  static const char *MAXTIME_PROP;
  static const char *ACTIVE_PROP;
  static const char *COARSEN_FACTOR;
  static const char *IMPOSED_STRAIN;

  enum StrainType
  {
    Free,
    ConstantRate,
    Path
  };

public: // Public Methods:
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
  FuncVector makeStrainFuncs_(const Vector &strainRate) const;
  FuncVector getStrainFuncs_(const Properties &globdat) const;
  void sortBndNodes_();
  template <typename T>
  void sortBndFace_(T &bndFace, const idx_t &index);
  void createTractionMesh_();
  void coarsenMesh_(FlexVector &trFace);
  void augmentFext_(const Vector &fext);
  void augmentMatrix_(Ref<MatrixBuilder> mbuilder, const Vector &fint, const Vector &disp);
  void getTractionMeshNodes_(IdxVector &connect, const Vector &x, const idx_t &face);

private: // Private Members (internal use)
  Assignable<XNodeSet> nodes_;
  BoolVector active_;
  idx_t rank_; // number of dimensions

  Ref<XDofSpace> dofs_;
  Ref<Constraints> cons_;
  IdxVector U_doftypes_;      // displacement dof types
  IdxVector T_doftypes_;      // traction dof types
  Ref<BoundaryShape> bshape_; // boundary element
  int nIP_;                   // number of int. points of boundary element
  int nnod_;                  // number of nodes of boundary element
  int localrank_;             // local rank of boundary element

  IdxVector bndNodes_[6];   // boundary nodes [ xmin, xmax, ymin, ymax ]
  FlexVector trNodes_[3];   // traction mesh nodes [ xmin, ymin ]
  Tuple<idx_t, 3> masters_; // corner nodes [ cornerX, cornerY, cornerZ ]
  idx_t ifixed_;            // master corner node [ corner0 ]

  Vector strain_; // applied Strain

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