/*
 *
 * A model to impose periodic boundary conditions on
 * the outer boundary of unit cells.
 *
 * Frans van der Meer, November 2014
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

#include "LagrangePeriodicModel.h"
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

//=========================================================
//    class LagrangePeriodicModel
//=========================================================
//-----------------------------------------------------------------------
//   static constants
//-----------------------------------------------------------------------

const char *LagrangePeriodicModel::STRAINRATE_PROP = "strainRate";
const char *LagrangePeriodicModel::STRAINPATH_PROP = "strainPath";
const char *LagrangePeriodicModel::MAXTIME_PROP = "maxTime";
const char *LagrangePeriodicModel::ACTIVE_PROP = "active";
const char *LagrangePeriodicModel::DUPEDNODES_PROP = "duplicatedNodes";

//-----------------------------------------------------------------------
//   constructor
//-----------------------------------------------------------------------

LagrangePeriodicModel::LagrangePeriodicModel

    (const String &name,
     const Properties &conf,
     const Properties &props,
     const Properties &globdat) : PeriodicBCModel::PeriodicBCModel(name, conf, props, globdat)
{
    System::warn() << "More creation stuff\n";
}

LagrangePeriodicModel::~LagrangePeriodicModel()
{
}

//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------

bool LagrangePeriodicModel::takeAction

    (const String &action,
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

        // apply strain (on corner nodes)

        if (strainType_ != Free)
        {
            PeriodicBCModel::applyStrain_(imposedStrain_);
        }

        // apply periodic conditions

        PeriodicBCModel::setPeriodicCons_();

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

void LagrangePeriodicModel::init_

    (const Properties &globdat)

{
    dofs_ = XDofSpace::get(nodes_.getData(), globdat);
    cons_ = Constraints::get(dofs_, globdat);

    dofTypes_.resize(rank_);
    dx_.resize(rank_);

    dofTypes_[0] = dofs_->addType("dx");
    dofTypes_[1] = dofs_->addType("dy");
    if (rank_ > 2)
    {
        dofTypes_[2] = dofs_->addType("dz");
    }

    // get nodes on boundaries

    for (idx_t i = 0; i < 2 * rank_; ++i)
    {
        NodeGroup edge = NodeGroup::get(PBCGroupInputModule::EDGES[i], nodes_, globdat, getContext());
        bndNodes_[i].ref(edge.getIndices());
    }
    for (idx_t i = 0; i < rank_; ++i)
    {
        if (bndNodes_[i * 2].size() != bndNodes_[i * 2 + 1].size())
        {
            throw Error(JEM_FUNC, String(i) +
                                      " opposite edges do not have the same number of nodes!!!");
        }
    }

    // get specimen dimensions
    // NB: assuming rectangular or hexahedral orientation
    //     this assumption is not checked in the code!

    Matrix coords(nodes_.toMatrix());
    Vector box(rank_ * 2);
    for (idx_t ix = 0; ix < rank_; ++ix)
    {
        box[ix * 2] = coords(ix, bndNodes_[ix * 2][0]);
        box[ix * 2 + 1] = coords(ix, bndNodes_[ix * 2 + 1][0]);

        dx_[ix] = box[ix * 2 + 1] - box[ix * 2];
    }

    // find which nodes are in the corners

    ifixed_ = NodeGroup::get(PBCGroupInputModule::CORNERS[0], nodes_, globdat, getContext()).getIndices()[0];

    for (idx_t i = 0; i < rank_; ++i)
    {
        masters_[i] = NodeGroup::get(PBCGroupInputModule::CORNERS[i + 1], nodes_, globdat, getContext()).getIndices()[0];
    }
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newLagrangePeriodicModel
//-----------------------------------------------------------------------

static Ref<Model> newLagrangePeriodicModel

    (const String &name,
     const Properties &conf,
     const Properties &props,
     const Properties &globdat)

{
    return newInstance<LagrangePeriodicModel>(name, conf, props, globdat);
}

//-----------------------------------------------------------------------
//   declareLagrangePeriodicModel
//-----------------------------------------------------------------------

void declareLagrangePeriodicModel()

{
    using jive::model::ModelFactory;

    ModelFactory::declare("LagPeriodicBC", &newLagrangePeriodicModel);
}
