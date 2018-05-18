#ifndef LAG_PERIODICBC_MODEL_H
#define LAG_PERIODICBC_MODEL_H

//-----------------------------------------------------------------------
//   class LagrangePeriodicModel
//-----------------------------------------------------------------------

#include "PeriodicBCModel.h"

#include <jem/util/Flex.h>
#include <jive/fem/XNodeSet.h>
#include <jive/algebra/MatrixBuilder.h>
#include <jive/geom/BoundaryLine.h>

using jem::util::Flex;
using jive::fem::XNodeSet;
using jive::algebra::MatrixBuilder;
using jive::geom::BoundaryLine2;
using jive::geom::BoundaryShape;

class LagrangePeriodicModel : public PeriodicBCModel
{
    public:
        typedef PeriodicBCModel Self;
        typedef Model Super;
        typedef Flex<idx_t> FlexVector;
        typedef Flex<idx_t>::Iterator Iter;
        typedef Array<Ref<Function>, 1> FuncVector;

        static const char *NODE_GROUPS[6];
        static const char *STRAINRATE_PROP;
        static const char *STRAINPATH_PROP;
        static const char *MAXTIME_PROP;
        static const char *ACTIVE_PROP;
        static const char *DUPEDNODES_PROP;
        static const char *IMPOSED_STRAIN;

        static const char *COARSEN_FACTOR;

        LagrangePeriodicModel

            (const String &name,
             const Properties &conf,
             const Properties &props,
             const Properties &globdat);

        virtual bool takeAction

            (const String &action,
             const Properties &params,
             const Properties &globdat);

    protected:
        virtual ~LagrangePeriodicModel();

        void init_(const Properties &globdat);

        void sortBndNodes_();
        template <typename T>
        void sortBndFace_(T &bndFace, const idx_t &index);
        void createTractionMesh_();
        void coarsenMesh_(FlexVector &trFace);
        void augmentFext_(const Vector &fext, const Vector &disp);
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
        double factor; // coarsening factor
};

#endif
