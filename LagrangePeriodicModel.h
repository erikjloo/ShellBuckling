#ifndef LAG_PERIODICBC_MODEL_H
#define LAG_PERIODICBC_MODEL_H

//-----------------------------------------------------------------------
//   class LagrangePeriodicModel
//-----------------------------------------------------------------------

#include "PeriodicBCModel.h"

class LagrangePeriodicModel : public PeriodicBCModel
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

    LagrangePeriodicModel

        (const String &name,
         const Properties &conf,
         const Properties &props,
         const Properties &globdat);

    // virtual void configure

    //     (const Properties &props,
    //      const Properties &globdat);

    // virtual void getConfig

    //     (const Properties &conf,
    //      const Properties &globdat) const;

    virtual bool takeAction

        (const String &action,
         const Properties &params,
         const Properties &globdat);

  protected:
    virtual ~LagrangePeriodicModel();

  protected:
    void init_

        (const Properties &globdat);

//     void setPeriodicCons_() const;

//     void advance_() const;

//     void checkCommit_

//         (const Properties &params,
//          const Properties &globdat) const;

//     void fixCorner_() const;

//     void applyStrain_

//         (const Vector &strain) const;

//     FuncVector makeStrainFuncs_

//         (const Vector &strainRate) const;

//     FuncVector getStrainFuncs_

//         (const Properties &globdat) const;

//   private:
//     IdxVector dofTypes_;

//     Assignable<NodeSet> nodes_;
//     BoolVector active_;

//     Ref<XDofSpace> dofs_;
//     Ref<Constraints> cons_;

//     IdxVector bndNodes_[6];   // indices of boundary nodes
//     Tuple<idx_t, 3> masters_; // master corner nodes
//     idx_t ifixed_;

//     Vector imposedStrain_; //total applied strain

//     double time_;
//     double stepSize_;
//     double maxTime_;

//     Vector dx_;
//     idx_t rank_;

//     String strainFile_;
//     StrainType strainType_;
//     FuncVector strainFunc_;

//     String dupedNodeGroup_;
};

#endif
