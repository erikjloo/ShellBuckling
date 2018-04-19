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
#include <jem/base/limits.h>
#include <jem/base/System.h>
#include <jem/base/Error.h>
#include <jem/numeric/func/UserFunc.h>
#include <jem/util/Properties.h>

#include <jive/model/Actions.h>
#include <jive/util/utilities.h>
#include <jive/util/FuncUtils.h>
#include <jive/util/Printer.h>

#include "WeakPeriodicBCModel.h"
#include "models.h"
#include "SolverNames.h"
#include "utilities.h"
#include "voigtUtilities.h"
#include "PBCGroupInputModule.h"

using jem::Error;
using jem::io::endl;
using jem::numeric::UserFunc;
using jive::fem::NodeGroup;
using jive::util::FuncUtils;

//=======================================================================
//   class WeakPeriodicBCModel
//=======================================================================

//-----------------------------------------------------------------------
//   static constants
//-----------------------------------------------------------------------

const char *WeakPeriodicBCModel::STRAINRATE_PROP = "strainRate";
const char *WeakPeriodicBCModel::STRAINPATH_PROP = "strainPath";
const char *WeakPeriodicBCModel::MAXTIME_PROP = "maxTime";
const char *WeakPeriodicBCModel::ACTIVE_PROP = "active";

//-----------------------------------------------------------------------
//   constructor
//-----------------------------------------------------------------------

WeakPeriodicBCModel::WeakPeriodicBCModel

    (const String &name,
     const Properties &conf,
     const Properties &props,
     const Properties &globdat) : Super(name)

{
  maxTime_ = jem::maxOf(maxTime_);

  Properties myProps = props.getProps(myName_);
  Properties myConf = conf.makeProps(myName_);

  // NodeSet and rank
  nodes_ = XNodeSet::find(globdat);
  rank_ = nodes_.rank();

  // Corner Nodes:
  masters_ = -1; // Initialize cornerx, cornery, cornerz
  ifixed_ = -1;  // Initialize corner 0

  // Variables required by applyStrain_()
  idx_t strCount = STRAIN_COUNTS[rank_];

  imposedStrain_.resize(strCount);
  imposedStrain_ = 0.;

  strainFunc_.resize(strCount);
  strainType_ = Free;

  Vector dStrain;
  myProps.find(dStrain, STRAINRATE_PROP);
  myConf.set(STRAINRATE_PROP, dStrain);

  bool strainPath = false;
  myProps.find(strainPath, STRAINPATH_PROP);
  myConf.set(STRAINPATH_PROP, strainPath);

  JEM_ASSERT(dStrain.size() == 0 || dStrain.size() == strCount);

  if (dStrain.size() == strCount)
  {
    strainType_ = ConstantRate;
    // oldStrain_     . resize ( strCount );
    // oldStrain_     = 0.;

    strainFunc_ = makeStrainFuncs_(dStrain);

    if (strainPath)
    {
      System::warn() << "strainPath input ignored, because strainRate "
                     << "has already been defined." << endl;
    }
  }
  else
  {
    if (strainPath)
    {
      strainType_ = Path;
      strainFunc_ = getStrainFuncs_(globdat);

      myProps.get(maxTime_, MAXTIME_PROP);
      myConf.set(MAXTIME_PROP, maxTime_);
    }
  }

  // Other variables
  if (myProps.find(active_, ACTIVE_PROP))
  {
    JEM_PRECHECK(active_.size() == rank_);
  }
  else
  {
    active_.resize(rank_);
    active_ = true;
  }

  time_ = 0.;
  stepSize_ = 1.;
}

//-----------------------------------------------------------------------
//   Destructor
//-----------------------------------------------------------------------

WeakPeriodicBCModel::~WeakPeriodicBCModel() {}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------

void WeakPeriodicBCModel::configure(const Properties &props,
                                    const Properties &globdat)
{
}

//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------

void WeakPeriodicBCModel::getConfig(const Properties &conf,
                                    const Properties &globdat) const
{
}

//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------

bool WeakPeriodicBCModel::takeAction(const String &action,
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
    advance_();
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
    fixCorner_();
    if (strainType_ != Free)
      applyStrain_(imposedStrain_);
    setPeriodicCons_();
    return true;
  }
  if (action == SolverNames::CHECK_COMMIT)
  {
    checkCommit_(params, globdat);
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

void WeakPeriodicBCModel::init_(const Properties &globdat)
{
  // Get dofs and cons
  dofs_ = XDofSpace::get(nodes_.getData(), globdat);
  cons_ = Constraints::get(dofs_, globdat);

  // Add displacement dof types
  U_doftypes_.resize(rank_);
  U_doftypes_[0] = dofs_->addType("dx");
  U_doftypes_[1] = dofs_->addType("dy");
  if (rank_ > 2)
    U_doftypes_[2] = dofs_->addType("dz");

  // Add traction dof types

  // Get boundary nodes
  for (idx_t face = 0; face < 2 * rank_; ++face)
  {
    NodeGroup edge = NodeGroup::get(PBCGroupInputModule::EDGES[face], nodes_, globdat, getContext());
    bndNodes_[face].ref(edge.getIndices());
  }
  System::warn() << "================================================\n";
  System::warn() << " bndNodes = [xmin, xmax, ymin, ymax] = \n";
  sortBndNodes_();

  // Get specimen dimensions
  dx_.resize(rank_);
  dx0_.resize(rank_);
  box_.resize(rank_ * 2);
  Matrix coords(nodes_.toMatrix());
  for (idx_t ix = 0; ix < rank_; ++ix)
  {
    box_[ix * 2] = coords(ix, bndNodes_[ix * 2][0]);
    box_[ix * 2 + 1] = coords(ix, bndNodes_[ix * 2 + 1][0]);
    dx_[ix] = box_[ix * 2 + 1] - box_[ix * 2];
  }
  System::warn() << "================================================\n";
  System::warn() << "Box dimensions = " << box_ << "\n";
  System::warn() << "dx_ = " << dx_ << "\n";

  // Find corner nodes
  ifixed_ = NodeGroup::get(PBCGroupInputModule::CORNERS[0], nodes_, globdat, getContext()).getIndices()[0];
  System::warn() << "================================================\n";
  System::warn() << "ifixed = Corner0 = " << ifixed_ << "\n";
  for (idx_t ix = 0; ix < rank_; ++ix)
    masters_[ix] = NodeGroup::get(PBCGroupInputModule::CORNERS[ix + 1], nodes_, globdat, getContext()).getIndices()[0];
  System::warn() << "Masters = [CornerX, CornerY] = " << masters_ << "\n";

  // Create traction mesh
  System::warn() << "================================================\n";
  System::warn() << " trNodes = [xmin, ymin] = \n";
  createTractionMesh_();

  // Create boundary element
  String bscheme = "Gauss1";
  bshape_ = BoundaryLine2::getShape("line", bscheme);
  nIP_ = bshape_->integrationPointCount();
  nnod_ = bshape_->nodeCount();
  localrank_ = bshape_->localRank();
}

//-----------------------------------------------------------------------
//   advance_
//-----------------------------------------------------------------------

void WeakPeriodicBCModel::advance_() const
{
  for (idx_t i = 0; i < imposedStrain_.size(); ++i)
  {
    if (strainFunc_[i] != NIL)
    {
      imposedStrain_[i] = strainFunc_[i]->eval(time_);
    }
  }
}

//-----------------------------------------------------------------------
//   checkCommit_
//-----------------------------------------------------------------------

void WeakPeriodicBCModel::checkCommit_(const Properties &params,
                                       const Properties &globdat) const
{
  // terminate the computation if displacement exceeds maximum.
  if (time_ > maxTime_)
  {
    System::out() << myName_ << " says: TERMINATE because "
                  << " time > maxTime." << endl;
    params.set(SolverNames::TERMINATE, "sure");
  }
}

//-----------------------------------------------------------------------
//   fixCorner_
//-----------------------------------------------------------------------

void WeakPeriodicBCModel::fixCorner_() const
{
  // apply zero displacement on first corner
  for (idx_t i = 0; i < rank_; ++i)
  {
    idx_t idof = dofs_->getDofIndex(ifixed_, U_doftypes_[i]);
    cons_->addConstraint(idof);
  }
}

//-----------------------------------------------------------------------
//   applyStrain_
//-----------------------------------------------------------------------

void WeakPeriodicBCModel::applyStrain_(const Vector &strain) const
{
  // prescribe displacement on the master corner nodes
  Matrix eps(rank_, rank_);
  voigtUtilities::voigt2TensorStrain(eps, strain);

  for (idx_t i = 0; i < rank_; ++i)
  {
    for (idx_t j = 0; j < rank_; ++j)
    {
      idx_t ivoigt = voigtUtilities::voigtIndex(i, j, rank_);
      if (strainFunc_[ivoigt] != NIL)
      {
        idx_t idof = dofs_->getDofIndex(masters_[i], U_doftypes_[j]);
        cons_->addConstraint(idof, dx_[i] * eps(i, j));
      }
    }
  }
}

//-----------------------------------------------------------------------
//   setPeriodicCons_
//-----------------------------------------------------------------------

void WeakPeriodicBCModel::setPeriodicCons_() const

{
  // loop over right and top edge: connect to opposite edges
  // relative to corner movement

  IdxVector idofs(2);
  Vector coefs(2);
  coefs = 1.;

  // loop over face pairs (ix)
  for (idx_t ix = 0; ix < rank_; ix++)
  {
    if (active_[ix])
    {
      IdxVector inodes(bndNodes_[ix * 2]);     // bndFace_min
      IdxVector jnodes(bndNodes_[ix * 2 + 1]); // bndFace_max

      // Loop over indices of inodes (bndFace_min)
      for (idx_t in = 0; in < inodes.size(); in++)
      {
        if (jnodes[in] != masters_[ix])
        {
          // Loop over dof types
          for (idx_t jx = 0; jx < rank_; ++jx)
          {
            // Master dofs: opposite node and relevant corner node
            idofs[0] = dofs_->getDofIndex(inodes[in], U_doftypes_[jx]);
            idofs[1] = dofs_->getDofIndex(masters_[ix], U_doftypes_[jx]);

            // Dofs of dependent (slave) node jd
            idx_t jdof = dofs_->getDofIndex(jnodes[in], U_doftypes_[jx]);

            // Add constraint jdof = idof[0]*coef[0] + idof[1]*coef[1]
            cons_->addConstraint(jdof, idofs, coefs);
          }
        }
      }
    }
  }
  //cons_->printTo(jive::util::Printer::get());
}

//-----------------------------------------------------------------------
//   makeStrainFuncs_
//-----------------------------------------------------------------------

Array<Ref<Function>, 1> WeakPeriodicBCModel::makeStrainFuncs_

    (const Vector &strainRate) const

{
  // make Function objects that define strain components as linear
  // function of time

  FuncVector ret(strainRate.size());

  for (idx_t i = 0; i < strainRate.size(); ++i)
  {
    String args = "t";
    String expr = String(strainRate[i]) + String("*t");
    ret[i] = newInstance<UserFunc>(args, expr);
  }

  return ret;
}

//-----------------------------------------------------------------------
//   getStrainFuncs_
//-----------------------------------------------------------------------

Array<Ref<Function>, 1> WeakPeriodicBCModel::getStrainFuncs_

    (const Properties &globdat) const

{
  System::out() << myName_ << " gets strain path from globdat\n"
                << "supposing UserFuncs have been defined elsewhere." << endl;

  FuncVector ret(strainFunc_.size());

  String context = "getStrainFuncs_(globdat)";

  // normal strain components are optional to allow for uniaxial load
  // shear strain components are required to exclude rigid body modes

  if (ret.size() == 3)
  {
    ret[0] = FuncUtils::findFunc("strain_xx", globdat);
    ret[1] = FuncUtils::findFunc("strain_yy", globdat);
    ret[2] = FuncUtils::getFunc("strain_xy", globdat, context);
  }
  else
  {
    ret[0] = FuncUtils::findFunc("strain_xx", globdat);
    ret[1] = FuncUtils::findFunc("strain_yy", globdat);
    ret[2] = FuncUtils::findFunc("strain_zz", globdat);
    ret[3] = FuncUtils::getFunc("strain_xy", globdat, context);
    ret[4] = FuncUtils::getFunc("strain_yz", globdat, context);
    ret[5] = FuncUtils::getFunc("strain_zx", globdat, context);
  }

  return ret;
}

//-----------------------------------------------------------------------
//   sortBndNodes_
//-----------------------------------------------------------------------

void WeakPeriodicBCModel::sortBndNodes_()
{
  // Loop over faces of bndNodes_
  for (idx_t face = 0; face < 2 * rank_; ++face)
  {
    // Map face onto ix
    int ix = face / 2;
    // Get correct index
    idx_t index = ((ix == 0) ? 1 : 0);
    // Perform bubblesort on bndNodes_[face]
    sortBndFace_(bndNodes_[face], index);
    // Print to verify
    System::warn() << " bndNodes_[" << face << "] = " << bndNodes_[face] << '\n';
  }
}

//-----------------------------------------------------------------------
//   sortbndFace_
//-----------------------------------------------------------------------

template <typename T>
void WeakPeriodicBCModel::sortBndFace_(T &bndFace, const idx_t &index)
{
  Vector c0(3); // coordinate vector of node "0"
  Vector c1(3); // coordinate vector of node "1"

  // Bubblesort algorithm
  for (idx_t in = 0; in < bndFace.size(); ++in)
  {
    for (idx_t jn = 0; jn < bndFace.size() - 1 - in; ++jn)
    {
      // Get nodal coordinatess
      nodes_.getNodeCoords(c0, bndFace[jn]);
      nodes_.getNodeCoords(c1, bndFace[jn + 1]);
      if (c0[index] > c1[index])
      {
        // Swap indices without creating a new temp variable
        bndFace[jn + 1] = bndFace[jn] + bndFace[jn + 1]; // b = a+b
        bndFace[jn] = bndFace[jn + 1] - bndFace[jn];     // a = (a+b)-a = b
        bndFace[jn + 1] = bndFace[jn + 1] - bndFace[jn]; // b = (a+b)-b = a
      }
    }
  }
}

//-----------------------------------------------------------------------
//   createTractionMesh_
//-----------------------------------------------------------------------

void WeakPeriodicBCModel::createTractionMesh_()
{
  //---------------------------------------------------------------------------
  // Part 1: Find the smallest element dimensions along the x and y coordinates
  //---------------------------------------------------------------------------

  Vector c0(3);       // coordinate vector of node "0"
  Vector c1(3);       // coordinate vector of node "1"
  double dx;          // dimensions of current element
  dx0_ = dx_.clone(); // smallest element dimensions

  // Loop over faces of bndNodes_
  for (idx_t face = 0; face < 2 * rank_; ++face)
  {
    // Map face onto ix
    int ix = face / 2;

    // Get correct index
    idx_t index = ((ix == 0) ? 1 : 0);

    // Assign bndNodes_[face] to bndFace
    IdxVector bndFace(bndNodes_[face]);

    // Loop over indices of bndFace
    for (idx_t in = 0; in < bndFace.size() - 1; ++in)
    {
      // get nodal coordinates
      nodes_.getNodeCoords(c0, bndFace[in]);
      nodes_.getNodeCoords(c1, bndFace[in + 1]);

      // Calculate dx and compare to dx0
      dx = c1[index] - c0[index];
      dx0_[index] = ((dx < dx0_[index]) ? dx : dx0_[index]);
    }
  }

  //---------------------------------------------------------------------------
  // Part 2 : Map all nodes onto xmin and ymin, reorder them and coarsen the mesh
  //---------------------------------------------------------------------------

  FlexVector trFace; // space for traction face
  Vector coords(rank_); // coordinate vector

  // Loop over faces of trNodes_
  for (idx_t ix = 0; ix < rank_; ++ix)
  {
    IdxVector inodes(bndNodes_[2 * ix]);     // bndFace_min
    IdxVector jnodes(bndNodes_[2 * ix + 1]); // bndFace_max
    trFace.resize(inodes.size() + jnodes.size());

    // Loop over indices of inodes
    idx_t jn = 0;
    for (idx_t in = 0; in < inodes.size(); ++in)
    {
      nodes_.getNodeCoords(coords, inodes[in]);
      trFace[jn++] = nodes_.addNode(coords);
    }

    // Loop over indices of jnodes
    for (idx_t in = 0; in < jnodes.size(); ++in)
    {
      nodes_.getNodeCoords(coords, jnodes[in]);
      coords[ix] = box_[2*ix];
      trFace[jn++] = nodes_.addNode(coords);
    }

    // Get correct index
    idx_t index = ((ix == 0) ? 1 : 0);

    // Sorting trFace (works only for 2D)
    sortBndFace_(trFace, index);

    // Coarsen the mesh (works only for 2D)
    coarsenMesh_(trFace, index);

    // Assign trFace to trNodes_[ix]
    trNodes_[ix] = trFace;

    // Print to verify
    System::warn() << "trNodes_[" << ix << "] = " << trNodes_[ix] << "\n";
  }

  // Add dofs to traction mesh
  for (idx_t ix = 0; ix < rank_; ++ix)
  {
    FlexVector trFace(trNodes_[ix]);
    for (idx_t in = 0; in < trFace.size(); ++in)
    {
      //dofs_->addDofs(trFace[in], T_doftypes_);
    }
  }
}

//-----------------------------------------------------------------------
//   coarsenMesh_
//-----------------------------------------------------------------------

void WeakPeriodicBCModel::coarsenMesh_(FlexVector &trFace, const idx_t &index)
{
  factor = 0.5;
  Vector dx(dx0_);
  dx[0] /= factor;
  dx[1] /= factor;

  Vector c0(3); // coordinate vector of node "0"
  Vector c1(3); // coordinate vector of node "1"
  Vector cn(3); // coordinate vector of node "n"
  nodes_.getNodeCoords(cn,trFace.back());

  // // Loop over indices of trFace
  // for (idx_t in = 0; in < trFace.size(); ++in)
  // {

  //   // Get nodal coordinates
  //   nodes_.getNodeCoords(c0, trFace[in]);
  //   nodes_.getNodeCoords(c1, trFace[in + 1]);

  //   // Delete indices
  //   while (c1[index] - c0[index] < dx[index])
  //   {
  //     // trFace.erase(in + 1);

  //     // Delete the current node

  //     // Assign c1 to the next node in line
  //     nodes_.getNodeCoords(c1, trFace[in + 1]);
  //   }
  //   if (cn[index] - c1[index] < dx[index])
  //   {
  //     // delete all nodes up to but not including the last one
  //   }
  // }
}

//-----------------------------------------------------------------------
//   getTractionMeshNodes_
//-----------------------------------------------------------------------

void WeakPeriodicBCModel::getTractionMeshNodes_(IdxVector &connect,
                                                const Vector &x,
                                                const idx_t &face)
{
  Matrix coords(rank_, nnod_);

  // Map face onto ix
  int ix = face / 2;

  // Implementation for two dimensions
  if (rank_ == 2)
  {
    // Assign trNodes_[ix] to trFace
    FlexVector trFace(trNodes_[ix]);

    // Loop over indices of trFace
    for (idx_t in = 0; in < trFace.size() - 1; ++in)
    {
      // get coords of nodes in trFace[in:in+2]
      connect[0] = trFace[in];
      connect[1] = trFace[in + 1];
      nodes_.getSomeCoords(coords, connect);

      // Get correct index
      idx_t index = ((ix == 0) ? 1 : 0);

      // Check if c0[index] < x[index] < c1[index]
      if ((coords(index, 0) < x[index]) && (x[index] < coords(index, 1)))
      {
        System::warn() << coords(index, 0) << " < " << x[index] << " < " << coords(index, 1) << "\n";
        break;
      }
    }
  }
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

void WeakPeriodicBCModel::augmentMatrix_(Ref<MatrixBuilder> mbuilder,
                                         const Vector &force,
                                         const Vector &disp)
{
  System::warn() << "===========================================\n";
  System::warn() << "Augmenting Matrix !!\n";

  // Variables related to both U mesh and T mesh:
  IdxVector connect(nnod_);    // node indices; connect = [ node1 node2 ]
  Matrix coords(rank_, nnod_); // node coordinates; coords  = [ x[1] x[2] ], x = [x y z]'

  // Variables related to element on U mesh:
  IdxVector idofs(nnod_ * rank_); // Dof indices related to U
  Vector w(nIP_);                 // Integration weights [ jw[1] jw[2] ], jw[ip] =j[ip]*w[ip]
  Matrix N(nnod_, nIP_);          //  shape functions of U mesh: N = [ n[1] n[2] ], n[ip] = [n1 n2]'
  Matrix X(rank_, nIP_);          // global coordinates of IP; X = [ x[1] x[2] ], x[i] = [x y z]'

  // Variables related to element on T mesh:
  IdxVector jdofs(nnod_ * rank_); // Dof indices related to T
  Vector u(localrank_);           // local coordinates of given X[ip]; u = xi
  Matrix H(nnod_, nIP_);          // shape functions of T; H = [ h[1] h[2] ], h = [h1 h2]'

  // Matrix to be assembled: K[idofs,jdofs] += w[ip]*N[ip]*H[ip]
  Matrix Ke(nnod_ * rank_, nnod_ * rank_);

  // Loop over faces of bndNodes_
  for (idx_t face = 0; face < 2 * rank_; ++face)
  {
    // Assign bndNodes_[face] to bndFace
    IdxVector bndFace(bndNodes_[face]);

    // Loop over indices of bndFace
    for (idx_t in = 0; in < bndFace.size() - 1; ++in)
    {
      connect[0] = bndFace[in];
      connect[1] = bndFace[in + 1];

      // Get idofs, w and N from displacement mesh
      dofs_->getDofIndices(idofs, connect, U_doftypes_);
      nodes_.getSomeCoords(coords, connect);
      bshape_->getIntegrationWeights(w, coords);
      N = bshape_->getShapeFunctions();
      bshape_->getGlobalIntegrationPoints(X, coords);

      // Get jdofs and H from traction mesh
      getTractionMeshNodes_(connect, X(ALL, 0), face);
      dofs_->getDofIndices(jdofs, connect, T_doftypes_);
      nodes_.getSomeCoords(coords, connect);

      Ke = 0.0;
      for (idx_t ip = 0; ip < nIP_; ip++)
      {
        bshape_->getLocalPoint(u, X(ALL, ip), 0.1, coords);
        bshape_->evalShapeFunctions(H[ip], u);

        // Ke += matmul(H[ip]*N[ip])*w[ip]
      }
      //H[ip]*N[ip]
      // w*H*N
      //mbuild
      // Assemble Ke
    }
  }
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newWeakPeriodicBCModel
//-----------------------------------------------------------------------

static Ref<Model> newWeakPeriodicBCModel

    (const String &name,
     const Properties &conf,
     const Properties &props,
     const Properties &globdat)

{
  return newInstance<WeakPeriodicBCModel>(name, conf, props, globdat);
}

//-----------------------------------------------------------------------
//   DeclareWeakPeriodicBCModel
//-----------------------------------------------------------------------

void declareWeakPeriodicBCModel()
{
  using jive::model::ModelFactory;

  ModelFactory::declare("WeakPeriodicBC", &newWeakPeriodicBCModel);
}