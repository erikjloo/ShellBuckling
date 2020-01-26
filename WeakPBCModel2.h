/*
 *
 * A model to impose weak periodic boundary 
 * conditions on the 2D unit cells. A different coarsening strategy is applied
 *
 * Erik Giesen Loo, April 2018
 * Frans van der Meer, Jan 2020
 *
 */

#ifndef WEAKPBC_MODEL_H_2
#define WEAKPBC_MODEL_H_2

//-----------------------------------------------------------------------
//   class WeakPBCModel2
//-----------------------------------------------------------------------

#include "PeriodicBCModel.h"

#include <jem/io/PrintWriter.h>
#include <jem/util/Flex.h>
#include <jive/fem/XNodeSet.h>
#include <jive/fem/XElementSet.h>
#include <jive/algebra/MatrixBuilder.h>
#include <jive/geom/BoundaryLine.h>

using jem::io::PrintWriter;
using jem::util::Flex;
using jive::algebra::MatrixBuilder;
using jive::fem::NodeGroup;
using jive::fem::XElementSet;
using jive::fem::XNodeSet;
using jive::geom::BoundaryLine2;
using jive::geom::BoundaryShape;
using jive::geom::Shape;

//=========================================================
//    class WeakPBCModel2
//=========================================================

class WeakPBCModel2 : public PeriodicBCModel
{
public:
    typedef WeakPBCModel2 Self;
    typedef PeriodicBCModel Super;
    typedef Flex<idx_t> FlexVector;
    typedef Flex<idx_t>::Iterator Iter;
    static const char *COARSEN_FACTOR;
    static const char *NUM_TNODE_PROP;
    static const double EPS;
    static const double PI;

    WeakPBCModel2

        (const String &name,
         const Properties &conf,
         const Properties &props,
         const Properties &globdat);

    virtual bool takeAction

        (const String &action,
         const Properties &params,
         const Properties &globdat);

protected:
    virtual ~WeakPBCModel2();

    void init_(const Properties &globdat);
    void sortBndNodes_();
    void findSmallestElement_();
    void clearTractionMesh_();
    void createTractionMesh_(); // discretization by coarsening factor
    void createTractionMesh2_(); // discretization by number of elements
    void sortBndFace_(IdxVector &bndFace, const idx_t &index);
    void sortTrFace_(FlexVector &trFace, const idx_t &index);
    void coarsenMesh_(FlexVector &trFace, const idx_t &index);
    void augmentMatrix_(Ref<MatrixBuilder> mbuilder, const Vector &fint, const Vector &disp);
    void getTractionMeshNodes_(IdxVector &connect, const Vector &x, const idx_t &face);
    // void createAlignedTractionMesh_();
    // void updateboundary_(const Properties &globdat);
    // void updateConnectivity_(const Properties &globdat);
    // void makeEdgeMesh_(XElementSet &edgeElems, const NodeGroup &edgeGroup, const idx_t ix);
    // void updateEdge_(const idx_t ix);

protected:

    idx_t numTNode_; // number of traction nodes
    double minDist_; // minimum relative size of integration element
    double angle_;   // alignment angle (clockwise)
    Vector box_;     // vector of specimen coordinates
    Vector dx0_;     // smallest element size per dimension
    double cf_;      // coarsening factor

    // Ref<BoundaryShape> ushape_; // for interpolation of displacements
    Ref<BoundaryShape> bshape_; // for integration of tractions
    // Ref<Shape> line2_;          // for integration

    IdxVector tdofTypes_;            // dof type indices for tractions
    FlexVector tnodeIndices_[3];     // array of flexvectors of traction mesh node incides
    Assignable<XNodeSet> nodes_;     // Nodeset for displacement (and tracion) mesh nodes
    Assignable<XNodeSet> tnodes_;    // NodeSet for traction mesh nodes
    Assignable<XElementSet> telems_; // ElementSet for traction mesh elements

    idx_t unodeCount_;   // number of nodes of ushape_
    idx_t tnodeCount_;   // number of nodes associated with bshape_
    idx_t udofCount_;    // numder of dofs for displacements
    idx_t tdofCount_;    // number of dofs associated with bshape_
    idx_t localrank_;    // number of local dimensions of bshape_
    idx_t ipCount_;      // number of int. points associated with bshape_

    // Tuple<Assignable<XElementSet>, 4> edgeElems_;
    // IdxVector edgeCount_;
    // idx_t inodeCount_; // number of nodes of inbshape_
    // idx_t irank_;      // dimenion of inbshape_
    // idx_t localRank_;  // local dimenion of ushape_
    // Vector offset_;    // offset in coordinates due to alignment
    // String ischemeString_;

    // Array<Ref<BndInfo>> bndInfo_;
};

#endif
