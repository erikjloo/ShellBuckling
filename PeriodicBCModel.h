#ifndef PERIODICBC_MODEL_H
#define PERIODICBC_MODEL_H

//-----------------------------------------------------------------------
//   class PeriodicBCModel
//-----------------------------------------------------------------------

#include <jive/model/Model.h>
#include <jive/model/ModelFactory.h>

#include <jive/fem/NodeGroup.h>
#include <jive/util/Assignable.h>
#include <jive/util/Constraints.h>
#include <jive/util/XDofSpace.h>

using namespace jem;

using jem::numeric::Function;
using jem::util::Properties;
using jive::Vector;
using jive::IdxVector;
using jive::BoolVector;
using jive::fem::NodeSet;
using jive::model::Model;
using jive::util::Assignable;
using jive::util::Constraints;
using jive::util::XDofSpace;

class PeriodicBCModel : public Model
{
  public:
    typedef PeriodicBCModel Self;
    typedef Model Super;
    typedef Array<Ref<Function>, 1> FuncVector;

    static const char *NODE_GROUPS[6];
    static const char *STRAINRATE_PROP;
    static const char *STRAINPATH_PROP;
    static const char *MAXTIME_PROP;
    static const char *ACTIVE_PROP;
    static const char *DUPEDNODES_PROP;
    static const char *IMPOSED_STRAIN;

    enum StrainType
    {
        Free,
        ConstantRate,
        Path
    };

    PeriodicBCModel

        (const String &name,
         const Properties &conf,
         const Properties &props,
         const Properties &globdat);

    virtual void configure

        (const Properties &props,
         const Properties &globdat);

    virtual void getConfig

        (const Properties &conf,
         const Properties &globdat) const;

    virtual bool takeAction

        (const String &action,
         const Properties &params,
         const Properties &globdat);

  protected:
    virtual ~PeriodicBCModel();

  protected:
    void init_

        (const Properties &globdat);

    void setPeriodicCons_() const;

    void advance_() const;

    void checkCommit_

        (const Properties &params,
         const Properties &globdat) const;

    void fixCorner_() const;

    void applyStrain_

        (const Vector &strain) const;

    FuncVector makeStrainFuncs_

        (const Vector &strainRate) const;

    FuncVector getStrainFuncs_

        (const Properties &globdat) const;

  protected:
    IdxVector dofTypes_;

    Assignable<NodeSet> nodes_;
    BoolVector active_;

    Ref<XDofSpace> dofs_;
    Ref<Constraints> cons_;

    IdxVector bndNodes_[6];   // indices of boundary nodes
    Tuple<idx_t, 3> masters_; // master corner nodes
    idx_t ifixed_;

    Vector imposedStrain_; //total applied strain

    double time_;
    double stepSize_;
    double maxTime_;

    Vector dx_;
    idx_t rank_;

    String strainFile_;
    StrainType strainType_;
    FuncVector strainFunc_;

    String dupedNodeGroup_;
};

#endif

