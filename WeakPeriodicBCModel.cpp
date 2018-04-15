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

  // Create traction mesh
  System::warn() << "================================================\n";
  System::warn() << " trNodes = [xmin, ymin] = \n";
  createTractionMesh_();

  // Get specimen dimensions
  dx_.resize(rank_);
  Matrix coords(nodes_.toMatrix());
  Vector box(rank_ * 2);
  // loop over face pairs (ix)
  for (idx_t ix = 0; ix < rank_; ++ix)
  {
    box[ix * 2] = coords(ix, bndNodes_[ix * 2][0]);
    box[ix * 2 + 1] = coords(ix, bndNodes_[ix * 2 + 1][0]);
    dx_[ix] = box[ix * 2 + 1] - box[ix * 2];
  }
  System::warn() << "================================================\n";
  System::warn() << "Box dimensions = " << box << "\n";

  // Get corner nodes
  ifixed_ = NodeGroup::get(PBCGroupInputModule::CORNERS[0], nodes_, globdat, getContext()).getIndices()[0];
  System::warn() << "================================================\n";
  System::warn() << "ifixed = Corner0 = " << ifixed_ << "\n";
  // loop over face pairs (ix)
  for (idx_t ix = 0; ix < rank_; ++ix)
    masters_[ix] = NodeGroup::get(PBCGroupInputModule::CORNERS[ix + 1], nodes_, globdat, getContext()).getIndices()[0];

  System::warn() << "===========================================\n";
  System::warn() << "Masters = [CornerX, CornerY] = " << masters_ << "\n";
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
      IdxVector inodes(bndNodes_[ix * 2]);     // bndFace_min (Gamma+)
      IdxVector jnodes(bndNodes_[ix * 2 + 1]); // bndFace_max

      // Loop over nodes indices of inodes (bndFace_min)
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
    // ix = 0 for xmin & xmax, 1 for ymin and ymax
    int ix = floor(face / 2);
    // index = 1 if ix = 0, else index = 0
    idx_t index = ((ix == 0)? 1 : 0);
    // Perform bubble sort on bndNodes_[face]
    sortBndFace_(bndNodes_[face], index);
    // Print to verify
    System::warn() << " bndNodes_[" << face << "] = " << bndNodes_[face] << '\n';
  }
}

//-----------------------------------------------------------------------
//   sortbndFace_
//-----------------------------------------------------------------------

void WeakPeriodicBCModel::sortBndFace_(IdxVector &bndFace, const idx_t &index)
{
  // Coordinate vectors of node "1" and node "0"
  Vector c1(3);
  Vector c0(3);

  // BubbleSort algorithm
  for (idx_t in = 0; in < bndFace.size(); ++in)
  {
    for (idx_t jn = 0; jn < bndFace.size() - 1 - in; ++jn)
    {
      nodes_.getNodeCoords(c1, bndFace[jn + 1]); // c1 = [x1, y1, z1]
      nodes_.getNodeCoords(c0, bndFace[jn]);     // c0 = [x0, y0, z0]
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
  IdxVector trFace;

  // Loop over face pairs (ix)
  for (idx_t ix = 0; ix < rank_; ix++)
  {
    IdxVector inodes(bndNodes_[ix * 2]);     // bndFace_min (Gamma+)
    IdxVector jnodes(bndNodes_[ix * 2 + 1]); // bndFace_max
    trFace.resize(inodes.size() + jnodes.size());

    // Loop over node indices in inodes
    idx_t jn = 0;
    Vector coords(rank_);
    for (idx_t in = 0; in < inodes.size(); ++in)
    {
      nodes_.getNodeCoords(coords, inodes[in]);
      trFace[jn++] = nodes_.addNode(coords);
    }

    // Loop over node indices in jnodes
    for (idx_t in = 0; in < jnodes.size(); ++in)
    {
      nodes_.getNodeCoords(coords, jnodes[in]);
      trFace[jn++] = nodes_.addNode(coords);
    }

    // Sorting trFace (works only for 2D)
    idx_t index = ((ix == 0) ? 1 : 0);
    sortBndFace_(trFace, index);

    // Assign trFace to trNodes_[ix]
    trNodes_[ix].ref(trFace.clone());

    // Print to verify
    System::warn() << "trNodes_[" << ix << "] = " << trNodes_[ix] << "\n";
  }

  // Add dofs to traction mesh
  for (idx_t ix = 0; ix < rank_; ++ix)
  {
    IdxVector trFace(trNodes_[ix]);
    for (idx_t in = 0; in < trFace.size(); ++in)
    {
      dofs_->addDofs(trFace[in], T_doftypes_);
    }
  }
}

//-----------------------------------------------------------------------
//   getTractionMeshNodes_
//-----------------------------------------------------------------------

void WeakPeriodicBCModel::getTractionMeshNodes_(IdxVector &connect,
                                                const Vector &x,
                                                const idx_t &face)
{
  Matrix coords(rank_, nnod_);

  // Map face index onto ix index
  int ix = floor(face / 2); 

  // Implementation for two dimensions
  if (rank_ == 2)
  {
    // Assign trNodes_[ix] to trFace
    IdxVector trFace(trNodes_[ix]);

    // Loop over node indices in trFace
    for (idx_t in = 0; in < trFace.size(); ++in)
    {
      // get coords of nodes in trFace[in:in+1]
      connect[0] = trFace[in];
      connect[1] = trFace[in + 1];
      nodes_.getSomeCoords(coords, connect);

      if (ix == 0)
      {
        if ((coords(1, 0) < x[1]) && (x[1] < coords(1, 1)))
          System::warn() << coords(1, 0) << " < " << x[1] << " < " << coords(1, 1) << "\n";
        break;
      }
      if (ix == 1)
      {
        if ((coords(0, 0) < x[0]) && (x[1] < coords(0, 1)))
          System::warn() << coords(0, 0) << " < " << x[0] << " < " << coords(0, 1) << "\n";
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

  // Loop over faces of bndNodes_
  for (idx_t face = 0; face < 2 * rank_; ++face)
  {
    // Assign bndNodes_[face] to bndFace
    IdxVector bndFace(bndNodes_[face]);

    // Loop over nodes indices in bndFace
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

      // Get jdofs and H from T mesh
      getTractionMeshNodes_(connect, X[0], face);
      dofs_->getDofIndices(jdofs, connect, T_doftypes_);
      nodes_.getSomeCoords(coords, connect);
      for (idx_t ip = 0; ip < nIP_; ip++)
      {
        bshape_->getLocalPoint(u, X[ip], 0.1, coords);
        bshape_->evalShapeFunctions(H[ip], u);
      }

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