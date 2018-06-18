/*
 *
 * A model to impose weak periodic boundary 
 * conditions on the 2D unit cells.
 *
 * Erik Giesen Loo, April 2018
 *
 */

#ifndef LAG_PERIODICBC_MODEL_H
#define LAG_PERIODICBC_MODEL_H

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
//    class LagrangePeriodicBCModel
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
    Assignable<XNodeSet> nodes_;

    Ref<BoundaryShape> bshape_; // boundary element
    int nIP_;                   // number of int. points of boundary element
    int nnod_;                  // number of nodes of boundary element
    int ndof_;                  // numder of dofs per boundary element
    int localrank_;             // local rank of boundary element

    FlexVector trNodes_[3]; // traction mesh nodes [ xmin, ymin ]

    Vector box_;   // specimen coordinates
    Vector dx0_;   // smallest element dimensions
    double factor_; // coarsening factor

    IdxVector trdofs_; // debugging vector
};

#endif
