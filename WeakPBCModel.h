/*
 *
 * A model to impose weak periodic boundary 
 * conditions on the 2D unit cells.
 *
 * Erik Giesen Loo, April 2018
 *
 */

#ifndef WEAKPBC_MODEL_H
#define WEAKPBC_MODEL_H

//-----------------------------------------------------------------------
//   class WeakPBCModel
//-----------------------------------------------------------------------

#include "PeriodicBCModel.h"

#include <jem/util/Flex.h>
#include <jive/fem/XNodeSet.h>
#include <jive/algebra/MatrixBuilder.h>
#include <jive/geom/BoundaryLine.h>

using jem::util::Flex;
using jive::algebra::MatrixBuilder;
using jive::fem::XNodeSet;
using jive::geom::BoundaryLine2;
using jive::geom::BoundaryShape;

//=========================================================
//    class WeakPBCModel
//=========================================================

class WeakPBCModel : public PeriodicBCModel
{
public:
  typedef Flex<idx_t> FlexVector;
  typedef Flex<idx_t>::Iterator Iter;
  static const char *COARSEN_FACTOR;

  WeakPBCModel

      (const String &name,
       const Properties &conf,
       const Properties &props,
       const Properties &globdat);

  virtual bool takeAction

      (const String &action,
       const Properties &params,
       const Properties &globdat);

protected:
  virtual ~WeakPBCModel();

  void init_(const Properties &globdat);
  void sortBndNodes_();
  template <typename T>
  void sortBndFace_(T &bndFace, const idx_t &index);
  void findSmallestElement_();
  void createTractionMesh_();
  void coarsenMesh_(FlexVector &trFace, const idx_t &index);
  void augmentMatrix_(Ref<MatrixBuilder> mbuilder, const Vector &fint, const Vector &disp);
  void getTractionMeshNodes_(IdxVector &connect, const Vector &x, const idx_t &face);

protected:
  Assignable<XNodeSet> nodes_; // XNodeSet to create the traction mesh nodes

  Ref<BoundaryShape> bshape_; // BoundaryShape for calculating shape functions

  int nIP_;       // number of int. points associated with bshape_
  int nnod_;      // number of nodes associated with bshape_
  int ndof_;      // numder of dofs associated with bshape_
  int localrank_; // number of local dimensions of bshape_

  FlexVector trNodes_[3]; // array of flexvectors of traction mesh nodes

  Vector box_;    // vector of specimen coordinates
  Vector dx0_;    // smallest element size per dimension
  double cf_; // coarsening factor
};

#endif
