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
  dx_.resize(rank_);

  // Add traction dof types


  // Get boundary nodes
  for (idx_t i = 0; i < 2 * rank_; ++i)
  {
    NodeGroup edge = NodeGroup::get(PBCGroupInputModule::EDGES[i], nodes_, globdat, getContext());
    bndNodes_[i].ref(edge.getIndices());
  }
  System::warn() << "===========================================\n";
  System::warn() << " bndNodes = [xmin, xmax, ymin, ymax] = \n";
  sortBndNodes_();

  // Create traction mesh
  System::warn() << "===========================================\n";
  System::warn() << " trNodes = [xmin, ymin] = \n";
  createTractionMesh_();

  // Get specimen dimensions
  Matrix coords(nodes_.toMatrix());
  Vector box(rank_ * 2);
  for (idx_t ix = 0; ix < rank_; ++ix)
  {
    box[ix * 2] = coords(ix, bndNodes_[ix * 2][0]);
    box[ix * 2 + 1] = coords(ix, bndNodes_[ix * 2 + 1][0]);
    dx_[ix] = box[ix * 2 + 1] - box[ix * 2];
  }
  System::warn() << "===========================================\n";
  System::warn() << "Box dimensions = " << box << "\n";

  // Get corner nodes
  ifixed_ = NodeGroup::get(PBCGroupInputModule::CORNERS[0], nodes_, globdat, getContext()).getIndices()[0];
  System::warn() << "===========================================\n";
  System::warn() << "ifixed = Corner0 = " << ifixed_ << "\n";

  for (idx_t i = 0; i < rank_; ++i)
    masters_[i] = NodeGroup::get(PBCGroupInputModule::CORNERS[i + 1], nodes_, globdat, getContext()).getIndices()[0];

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

  // loop over faces
  for (idx_t ix = 0; ix < rank_; ix++)
  {
    if (active_[ix])
    {
      IdxVector inodes(bndNodes_[ix * 2]);
      IdxVector jnodes(bndNodes_[ix * 2 + 1]);

      // loop over nodes of edges

      for (idx_t in = 0; in < inodes.size(); in++)
      {
        if (jnodes[in] != masters_[ix])
        {
          // System::out() << "constraining node " << jnodes[in] <<
          // " to  " << inodes[in] << " and " << masters_[ix] << endl;

          // loop over dof types

          for (idx_t jx = 0; jx < rank_; ++jx)
          {
            // master dofs: opposite node and relevant corner node

            idofs[0] = dofs_->getDofIndex(inodes[in], U_doftypes_[jx]);
            idofs[1] = dofs_->getDofIndex(masters_[ix], U_doftypes_[jx]);

            // dofs of dependent (slave) node jd

            idx_t jdof = dofs_->getDofIndex(jnodes[in], U_doftypes_[jx]);

            // add constraint
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
  idx_t index;
  // loop over faces of bndNodes_ (faces)
  for (idx_t face = 0; face < 2 * rank_; ++face)
  {
    if ((face == 0) || (face == 1))
      index = 1; // Index to compare y coordinates
    if ((face == 2) || (face == 3))
      index = 0; // Index to compare x coordinates
    // Perform bubble sort on bndNodes_[i]
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
  for (idx_t i = 0; i < bndFace.size(); ++i)
  {
    for (idx_t j = 0; j < bndFace.size() - 1 - i; ++j)
    {
      nodes_.getNodeCoords(c1, bndFace[j + 1]); // c1 = [x1, y1, z1]
      nodes_.getNodeCoords(c0, bndFace[j]);     // c0 = [x0, y0, z0]
      if (c0[index] > c1[index])
      {
        // Swap indices without creating a new temp variable
        bndFace[j + 1] = bndFace[j] + bndFace[j + 1]; // b = a+b
        bndFace[j] = bndFace[j + 1] - bndFace[j];     // a = (a+b)-a = b
        bndFace[j + 1] = bndFace[j + 1] - bndFace[j]; // b = (a+b)-b = a
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
  // loop over face pairs
  // Recall: bndNodes_ = [ xmin, xmax, ymin, ymax, zmin, zmax ]
  for (idx_t face = 0; face < rank_; ++face)
  {
    IdxVector bndFace_min(bndNodes_[face * 2]);
    IdxVector bndFace_max(bndNodes_[face * 2 + 1]);
    trFace.resize(bndFace_min.size() + bndFace_max.size());

    // loop over node indices in bndFace_min
    idx_t jnod = 0;
    Vector coords(rank_);
    for (idx_t inod = 0; inod < bndFace_min.size(); ++inod)
    {
      nodes_.getNodeCoords(coords, bndFace_min[inod]);
      trFace[jnod++] = nodes_.addNode(coords);
    }

    // loop over node indices in bndFace_max
    for (idx_t inod = 0; inod < bndFace_max.size(); ++inod)
    {
      nodes_.getNodeCoords(coords, bndFace_max[inod]);
      trFace[jnod++] = nodes_.addNode(coords);
    }

    // Sorting trFace
    idx_t index = ((face == 0) ? 1 : 0);
    sortBndFace_(trFace, index);
    // Print to verify
    // trNodes_[face] = trFace.clone();
    System::warn() << "trNodes_[" << face << "] = " << trFace << "\n";
  }

  System::warn() << "===========================================\n";
  System::warn() << "Adding dofs to traction mesh !!\n";
  // Add dofs to traction mesh
  for (idx_t face = 0; face < rank_; ++face)
  {
    IdxVector trFace(trNodes_[face]);
    for (idx_t inod = 0; inod < trFace.size(); ++inod)
    {
      dofs_->addDofs(trFace[inod], T_doftypes_);
    }
  }
}

//-----------------------------------------------------------------------
//   augmentMatrix_
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
//   augmentMatrix_
//-----------------------------------------------------------------------

void WeakPeriodicBCModel::augmentMatrix_(Ref<MatrixBuilder> mbuilder,
                                         const Vector &force,
                                         const Vector &disp)
{

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
  // loop over faces of bndNodes_ (faces)
  for (idx_t face = 0; face < 2 * rank_; ++face)
  {
    IdxVector bndFace(bndNodes_[face]);

    // loop over nodes indices in bndFace
    for (idx_t in = 0; in < bndFace.size() - 1; ++in)
    {
      connect[0] = bndFace[in];
      connect[1] = bndFace[in + 1];

      dofs_->getDofIndices(idofs, connect, U_doftypes_);

      nodes_.getSomeCoords(coords, connect);

      bshape_->getIntegrationWeights(w, coords);

      N = bshape_->getShapeFunctions();

      bshape_->getGlobalIntegrationPoints(X, coords);

      // getTractionMeshNodes_(connect, X);

      connect[0] = 3;
      connect[1] = 4;

      dofs_->getDofIndices(jdofs, connect, T_doftypes_);

      nodes_.getSomeCoords(coords, connect);

      for (idx_t ip = 0; ip < nIP_; ip++)
      {
        bshape_->getLocalPoint(u, X[ip], 0.1, coords);
        bshape_->evalShapeFunctions(H[ip], u);
      }
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
