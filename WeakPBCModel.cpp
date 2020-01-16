/*
 *
 * A model to impose weak periodic boundary
 * conditions on the 2D unit cells.
 *
 * Erik Giesen Loo, April 2018
 *
 */

#include <jem/base/array/utilities.h>
#include <jem/base/array/operators.h>
#include <jem/base/array/select.h>
#include <jem/base/limits.h>
#include <jem/base/System.h>
#include <jem/base/Error.h>
#include <jem/numeric/func/UserFunc.h>
#include <jem/numeric/algebra.h>
#include <jem/util/Properties.h>

#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
#include <jive/util/utilities.h>
#include <jive/util/FuncUtils.h>
#include <jive/util/Printer.h>

#include "WeakPBCModel.h"
#include "models.h"
#include "SolverNames.h"
#include "utilities.h"
#include "voigtUtilities.h"
#include "PBCGroupInputModule.h"

using jem::Error;
using jem::io::endl;
using jem::numeric::matmul;
using jem::numeric::UserFunc;
using jive::fem::NodeGroup;
using jive::model::StateVector;
using jive::util::FuncUtils;

//=========================================================
//    class WeakPBCModel
//=========================================================
//-----------------------------------------------------------------------
//   static constants
//-----------------------------------------------------------------------

const char *WeakPBCModel::COARSEN_FACTOR = "coarsenFactor";

//-----------------------------------------------------------------------
//   constructor
//-----------------------------------------------------------------------

WeakPBCModel::WeakPBCModel

    (const String &name,
     const Properties &conf,
     const Properties &props,
     const Properties &globdat) : PeriodicBCModel::PeriodicBCModel(name, conf, props, globdat)
{
  Properties myProps = props.getProps(myName_);
  Properties myConf = conf.makeProps(myName_);

  myProps.get(cf_, COARSEN_FACTOR);
  myConf.set(COARSEN_FACTOR, cf_);
  nodes_ = XNodeSet::find(globdat);
}

//-----------------------------------------------------------------------
//   Destructor
//-----------------------------------------------------------------------

WeakPBCModel::~WeakPBCModel()
{
}

//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------

bool WeakPBCModel::takeAction(const String &action,
                                       const Properties &params,
                                       const Properties &globdat)
{
  using jive::model::ActionParams;
  using jive::model::Actions;

  if (action == Actions::INIT)
  {
    init_(globdat);
    return true;
  }

  if (action == Actions::ADVANCE)
  {
    PeriodicBCModel::advance_();
    return true;
  }

  if (action == SolverNames::SET_STEP_SIZE)
  {
    double dt, dt0;
    params.get(dt, SolverNames::STEP_SIZE);
    params.get(dt0, SolverNames::STEP_SIZE_0);
    stepSize_ = dt / dt0;
    return true;
  }

  if (action == Actions::GET_CONSTRAINTS)
  {
    PeriodicBCModel::fixCorner_();
    if (strainType_ != Free)
      PeriodicBCModel::applyStrain_(imposedStrain_);
    return true;
  }

  if (action == Actions::GET_MATRIX0 || action == Actions::GET_INT_VECTOR)
  {
    Ref<MatrixBuilder> mbuilder;
    Vector solu;
    Vector fint;

    // Get the current displacements.

    StateVector::get(solu, dofs_, globdat);

    // Get the matrix builder and the internal force vector.

    params.get(fint, ActionParams::INT_VECTOR);
    params.find(mbuilder, ActionParams::MATRIX0);

    augmentMatrix_(mbuilder, fint, solu);
    return true;
  }

  if (action == SolverNames::CHECK_COMMIT)
  {
    PeriodicBCModel::checkCommit_(params, globdat);
    return true;
  }

  if (action == Actions::COMMIT)
  {
    time_ += stepSize_;
    return true;
  }

  return false;
}

//-----------------------------------------------------------------------
//   init_
//-----------------------------------------------------------------------

void WeakPBCModel::init_(const Properties &globdat)
{
  // Get dofs and cons
  dofs_ = XDofSpace::get(nodes_.getData(), globdat);
  cons_ = Constraints::get(dofs_, globdat);

  dofTypes_.resize(rank_);
  box_.resize(rank_ * 2);
  dx0_.resize(rank_);
  dx_.resize(rank_);

  // Add displacement dof types
  dofTypes_[0] = dofs_->addType("dx");
  dofTypes_[1] = dofs_->addType("dy");
  if (rank_ > 2)
    dofTypes_[2] = dofs_->addType("dz");

  // Get boundary nodes
  for (idx_t face = 0; face < 2 * rank_; ++face)
  {
    NodeGroup edge = NodeGroup::get(PBCGroupInputModule::EDGES[face], nodes_, globdat, getContext());
    bndNodes_[face].ref(edge.getIndices());
  }
  System::out() << "bndNodes = [xmin, xmax, ymin, ymax] = \n";
  sortBndNodes_();

  // Get specimen dimensions
  Matrix coords(nodes_.toMatrix());
  for (idx_t ix = 0; ix < rank_; ++ix)
  {
    box_[ix * 2] = coords(ix, bndNodes_[ix * 2][0]);
    box_[ix * 2 + 1] = coords(ix, bndNodes_[ix * 2 + 1][0]);
    dx_[ix] = box_[ix * 2 + 1] - box_[ix * 2];
  }
  System::out() << "box_ = " << box_ << "\n";
  System::out() << "dx_ = " << dx_ << "\n";

  // Find corner nodes
  ifixed_ = NodeGroup::get(PBCGroupInputModule::CORNERS[0], nodes_, globdat, getContext()).getIndices()[0];
  System::out() << "ifixed_ = Corner0 = " << ifixed_ << "\n";
  for (idx_t i = 0; i < rank_; ++i)
    masters_[i] = NodeGroup::get(PBCGroupInputModule::CORNERS[i + 1], nodes_, globdat, getContext()).getIndices()[0];
  System::out() << "masters_ = [CornerX, CornerY] = " << masters_ << "\n";

  // Create traction mesh
  System::out() << " trNodes_ = [xmin, ymin] = \n";
  findSmallestElement_();
  createTractionMesh_();

  // Create boundary element
  String bscheme = "Gauss1";
  bshape_ = BoundaryLine2::getShape("line", bscheme);
  nIP_ = bshape_->integrationPointCount();
  nnod_ = bshape_->nodeCount();
  ndof_ = nnod_ * rank_;
  localrank_ = bshape_->localRank();
}

//-----------------------------------------------------------------------
//   sortBndNodes_
//-----------------------------------------------------------------------

void WeakPBCModel::sortBndNodes_()
{
  // Loop over faces of bndNodes_
  for (idx_t face = 0; face < 2 * rank_; ++face)
  {
    // Map face onto ix
    int ix = face / 2;
    // Get correct index
    idx_t iy = (ix + 1) % 2;
    // Perform bubblesort on bndNodes_[face]
    sortBndFace_(bndNodes_[face], iy);
    // Print to verify
    System::out() << " bndNodes_[" << face << "] = " << bndNodes_[face] << '\n';
  }
}

//-----------------------------------------------------------------------
//   sortbndFace_
//-----------------------------------------------------------------------

template <typename T>
void WeakPBCModel::sortBndFace_(T &bndFace, const idx_t &iy)
{
  Vector c0(rank_); // coordinate vector of node "0"
  Vector c1(rank_); // coordinate vector of node "1"

  // Bubblesort algorithm
  for (idx_t in = 0; in < bndFace.size(); ++in)
  {
    for (idx_t jn = 0; jn < bndFace.size() - 1 - in; ++jn)
    {
      // Get nodal coordinatess
      nodes_.getNodeCoords(c0, bndFace[jn]);
      nodes_.getNodeCoords(c1, bndFace[jn + 1]);

      // Swap indices if necessary
      if (c0[iy] > c1[iy])
      {
        bndFace[jn + 1] = bndFace[jn] + bndFace[jn + 1]; // b = a+b
        bndFace[jn] = bndFace[jn + 1] - bndFace[jn];     // a = (a+b)-a = b
        bndFace[jn + 1] = bndFace[jn + 1] - bndFace[jn]; // b = (a+b)-b = a
      }
    }
  }
}

//-----------------------------------------------------------------------
//   findSmallestElement_
//-----------------------------------------------------------------------

void WeakPBCModel::findSmallestElement_()
{
  Vector c0(rank_);   // coordinate vector of node "0"
  Vector c1(rank_);   // coordinate vector of node "1"
  double dx;          // dimensions of current element
  dx0_ = dx_.clone(); // smallest element dimensions

  // Loop over faces of bndNodes_
  for (idx_t face = 0; face < 2 * rank_; ++face)
  {
    // Map face onto ix
    int ix = face / 2;

    // Get correct index
    idx_t iy = (ix + 1) % 2;

    // Loop over indices of bndFace
    for (idx_t in = 0; in < bndNodes_[face].size() - 1; ++in)
    {
      // Get nodal coordinates
      nodes_.getNodeCoords(c0, bndNodes_[face][in]);
      nodes_.getNodeCoords(c1, bndNodes_[face][in + 1]);

      // Calculate dx and compare to dx0
      dx = c1[iy] - c0[iy];
      dx0_[iy] = ((dx < dx0_[iy]) ? dx : dx0_[iy]);
    }
  }
}

//-----------------------------------------------------------------------
//   createTractionMesh_
//-----------------------------------------------------------------------

void WeakPBCModel::createTractionMesh_()
{
  Vector coords(rank_); // coordinate vector

  // Loop over faces of trNodes_
  for (idx_t ix = 0; ix < rank_; ++ix)
  {
    IdxVector inodes(bndNodes_[2 * ix]);     // bndFace_min
    IdxVector jnodes(bndNodes_[2 * ix + 1]); // bndFace_max
    trNodes_[ix].resize(inodes.size() + jnodes.size());

    // Loop over indices of inodes
    idx_t kn = 0;
    for (idx_t in = 0; in < inodes.size(); ++in)
    {
      nodes_.getNodeCoords(coords, inodes[in]);
      coords[ix] = box_[2 * ix + 1];
      trNodes_[ix][kn++] = nodes_.addNode(coords);
    }

    // Loop over indices of jnodes
    for (idx_t jn = 0; jn < jnodes.size(); ++jn)
    {
      nodes_.getNodeCoords(coords, jnodes[jn]);
      trNodes_[ix][kn++] = nodes_.addNode(coords);
    }

    // Get correct index
    idx_t iy = (ix + 1) % 2;

    // Sorting trFace (works only for 2D)
    sortBndFace_(trNodes_[ix], iy);

    // Coarsen the mesh (works only for 2D)
    coarsenMesh_(trNodes_[ix], iy);

    // if (ix == rank_-1)
    //   trNodes_[ix].back() = trNodes_[ix-1].back();

    // Add dofs to trFace
    for (idx_t in = 0; in < trNodes_[ix].size(); ++in)
      dofs_->addDofs(trNodes_[ix][in], dofTypes_);

    // Print to verify
    System::out() << "trNodes_[" << ix << "] = " << trNodes_[ix] << "\n";
  }
}

//-----------------------------------------------------------------------
//   coarsenMesh_
//-----------------------------------------------------------------------

void WeakPBCModel::coarsenMesh_(FlexVector &trFace, const idx_t &iy)
{
  using jem::numeric::norm2;

  Vector c0(rank_); // coordinate vector of node "0"
  Vector c1(rank_); // coordinate vector of node "1"
  Vector cn(rank_); // coordinate vector of node "n"
  nodes_.getNodeCoords(cn, trFace.back());
  
  // Minimum coarsening factor for Neumann boundary conditions
  double cf_min = (dx0_[0]+dx0_[1])/rank_/max(dx_[0],dx_[1]);
  cf_ = (cf_ < cf_min)? cf_min : cf_;
  // System::out() << "Minimum coarsening cf is " << cf_min << "\n";

  // Tolerance for removing nodes
  double dx = (dx0_[0] + dx0_[1]) / (2 * cf_);

  // Loop over indices of trFace
  int jn = 0;
  for (Iter in = trFace.begin(); in < trFace.end()-1; ++in)
  {
    // Get nodal coordinates
    nodes_.getNodeCoords(c0, trFace[jn]);
    nodes_.getNodeCoords(c1, trFace[jn + 1]);

    // Delete indices until c1 - c0 > dx
    while (norm2(c0 - c1) < min(dx, dx_[iy]) && jn < trFace.size())
    {
      // Delete current node index
      trFace.erase(in + 1);
      // Assign next nodal coordinates to c1
      nodes_.getNodeCoords(c1, trFace[jn + 1]);
    }

    // Check distance to last node
    if (norm2(cn - c1) < dx)
    {
      // Delete all indices up to but not including the last one
      trFace.erase(in + 1, trFace.end() - 1);
      // Break for loop of indices of trFace
      return;
    }
    jn += 1;
  }
}

//-----------------------------------------------------------------------
//   getTractionMeshNodes_
//-----------------------------------------------------------------------

void WeakPBCModel::getTractionMeshNodes_(IdxVector &connect,
                                        const Vector &x,
                                        const idx_t &face)
{
  Matrix coords(rank_, nnod_);

  // Map face onto ix
  int ix = face / 2;

  // Implementation for two dimensions
  if (rank_ == 2)
  {
    // Loop over indices of trFace
    for (idx_t in = 0; in < trNodes_[ix].size() - 1; ++in)
    {
      // get coords of nodes in trFace[in:in+2]
      connect[0] = trNodes_[ix][in];
      connect[1] = trNodes_[ix][in + 1];
      nodes_.getSomeCoords(coords, connect);

      // Get correct index
      idx_t iy = (ix + 1) % 2;

      // Check if c0[iy] < x[iy] < c1[iy]
      if ((coords(iy, 0) <= x[iy]) && (x[iy] <= coords(iy, 1)))
      {
        return;
      }
    }
  }

  throw jem::Error(JEM_FUNC, "nodes not found");
  // Implementation for three dimensions
  // if (rank_ == 3)
  // IdxVector trFace(trElems_[ix]);
  //
  // Loop over elements in trFace
  // for (idx_t ie = 0; ie < trFace.size(); ++ie)
  // {
  //    elems_.getNodeIndices(connect, ie);
  //    nodes_.getSomeCoords(coords, connect);
  //
  //    if (ix == 0)
  //    {
  //
  //    }
  //    if (ix == 1)
  //    {
  //
  //    }
  //    if (ix == 2)
  //    {
  //
  //    }
  // }
}

//-----------------------------------------------------------------------
//   augmentMatrix_
//-----------------------------------------------------------------------

void WeakPBCModel::augmentMatrix_(Ref<MatrixBuilder> mbuilder,
                                           const Vector &fint,
                                           const Vector &solu)
{
  // Variables related to both U mesh and T mesh:
  IdxVector connect(nnod_);    // node indices = [ node1 node2 ]
  Matrix coords(rank_, nnod_); // node coordinates = [ x[1] x[2] ], x = [x y z]'

  // Variables related to element on U mesh:
  IdxVector idofs(ndof_); // dof indices related to U
  Vector w(nIP_);         // Integration weights = [ jw[1] jw[2] ]
  Matrix n(nnod_, nIP_);  // shape functions of U mesh: N = [ n[1] n[2] ], n[ip] = [n1 n2]'
  Matrix N(rank_, ndof_); // N matrix of U mesh [n1 0 n2 0; 0 n1 0 n2]
  Matrix X(rank_, nIP_);  // global coordinates of IP; X = [ x[1] x[2] ], x[ip] = [x y z]'
  N = 0.0;

  // Variables related to element on T mesh:
  IdxVector jdofs(ndof_); // dof indices related to T
  Vector xi(localrank_);  // local coordinates of given X[ip]; u = xi
  Vector h(nnod_);        // shape functions of T; h = [h1 h2]'
  Matrix H(rank_, ndof_); // H matrix of T mesh [h1 0 h2 0; 0 h1 0 h2]
  H = 0.0;

  // Matrices and vector to be assembled:
  Vector tr(ndof_);         // Vector of tractions of each element
  Vector u(ndof_);          // Vector of displacements of each element
  Matrix Ke(ndof_, ndof_);  // Ke = w[ip]*N[ip]*H[ip]
  Matrix KeT(ndof_, ndof_); // KeT = w[ip]*H[ip]*N[ip]
  Vector fe(ndof_);         // fe = Ke*t or KeT*u

  // Loop over faces of bndNodes_
  for (idx_t face = 0; face < 2 * rank_; ++face)
  {
    // Loop over indices of bndFace
    for (idx_t in = 0; in < bndNodes_[face].size() - 1; ++in)
    {
      // Get idofs and w from U-mesh
      connect[0] = bndNodes_[face][in];
      connect[1] = bndNodes_[face][in + 1];
      nodes_.getSomeCoords(coords, connect);
      bshape_->getIntegrationWeights(w, coords);
      dofs_->getDofIndices(idofs, connect, dofTypes_);
      bshape_->getGlobalIntegrationPoints(X, coords);
      n = bshape_->getShapeFunctions();

      for (idx_t ip = 0; ip < nIP_; ip++)
      {
        // Get jdofs from T-mesh
        getTractionMeshNodes_(connect, X[ip], face);
        dofs_->getDofIndices(jdofs, connect, dofTypes_);
        nodes_.getSomeCoords(coords, connect);

        // Get N-matrix from U-mesh
        N(0, 0) = N(1, 1) = n(0, ip);
        N(0, 2) = N(1, 3) = n(1, ip);

        // Get H-matrix from T-mesh
        int ix = face / 2;
        X(ix,ip) = box_[2*ix+1];
        double eps = dx0_[ix]/100;
        bshape_->getLocalPoint(xi, X[ip], eps, coords);
        bshape_->evalShapeFunctions(h, xi);
        H(0, 0) = H(1, 1) = h[0];
        H(0, 2) = H(1, 3) = h[1];

        // Assemble Ke
        if (face == 0 || face == 2 || face == 4)
          Ke = +w[ip] * matmul(N.transpose(), H);
        else if (face == 1 || face == 3 || face == 5)
          Ke = -w[ip] * matmul(N.transpose(), H);

        // Add Ke and KeT to mbuilder
        KeT = Ke.transpose();
        if (mbuilder != NIL)
        {
          mbuilder->addBlock(idofs, jdofs, Ke);
          mbuilder->addBlock(jdofs, idofs, KeT);
        }

        // Assemble U-mesh fe
        tr = select(solu, jdofs);
        select(fint, idofs) += matmul(Ke, tr);

        // Assemble T-mesh fe
        u = select(solu, idofs);
        select(fint, jdofs) += matmul(KeT, u);
      }
    }
  }

  // Variables related to corner displacements
  Matrix Ht(ndof_, rank_);
  Vector u_corner(rank_);
  Vector u_fixed(rank_);
  IdxVector kdofs(rank_);
  dofs_->getDofIndices(kdofs, ifixed_, dofTypes_);
  u_fixed = solu[kdofs];

  // Loop over faces of trNodes_
  for (idx_t ix = 0; ix < rank_; ++ix)
  {
    // Obtain u_corner
    IdxVector idofs(rank_);
    dofs_->getDofIndices(idofs, masters_[ix], dofTypes_);
    u_corner = solu[idofs];

    // Loop over indices of trFace
    for (idx_t jn = 0; jn < trNodes_[ix].size() - 1; ++jn)
    {
      // Get jdofs, w and H from traction mesh
      connect[0] = trNodes_[ix][jn];
      connect[1] = trNodes_[ix][jn + 1];
      nodes_.getSomeCoords(coords, connect);
      bshape_->getIntegrationWeights(w, coords);
      dofs_->getDofIndices(jdofs, connect, dofTypes_);
      n = bshape_->getShapeFunctions();

      for (idx_t ip = 0; ip < nIP_; ip++)
      {
        // Assemble H matrix
        H(0, 0) = H(1, 1) = n(0, ip);
        H(0, 2) = H(1, 3) = n(1, ip);
        H *= w[ip];

        // Add H and Ht to mbuilder
        Ht = H.transpose();
        if (mbuilder != NIL)
        {
          mbuilder->addBlock(idofs, jdofs, H);
          mbuilder->addBlock(jdofs, idofs, Ht);
        }
        // Assemble U-mesh and T-mesh fe
        tr = select(solu, jdofs);
        select(fint, idofs) += matmul(H, tr);
        select(fint, jdofs) += matmul(Ht, u_corner);

        H = -H;
        Ht = -Ht;
        if (mbuilder != NIL)
        {
          mbuilder->addBlock(kdofs, jdofs, H);
          mbuilder->addBlock(jdofs, kdofs, Ht);
        }
        // Assemble U-mesh and T-mesh fe
        select(fint, kdofs) += matmul(H, tr);
        select(fint, jdofs) += matmul(Ht, u_fixed);
      }
    }
  }
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newWeakPBCModel
//-----------------------------------------------------------------------

static Ref<Model> newWeakPBCModel

    (const String &name,
     const Properties &conf,
     const Properties &props,
     const Properties &globdat)

{
  return newInstance<WeakPBCModel>(name, conf, props, globdat);
}

//-----------------------------------------------------------------------
//   declareWeakPBCModel
//-----------------------------------------------------------------------

void declareWeakPBCModel()

{
  using jive::model::ModelFactory;
  ModelFactory::declare("WeakPBC", &newWeakPBCModel);
}
