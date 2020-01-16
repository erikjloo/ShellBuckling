/*
 *
 * A model to impose weak periodic boundary
 * conditions on the 2D unit cells.
 *
 * Erik Giesen Loo, Jan 2020
 * Frans van der Meer, July 2019
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
// #include <jem/util/Event.h>

// #include <jive/fem/NodeGroup.h>
#include <jive/geom/Line.h>
// #include <jive/geom/ShapeFactory.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
// #include <jive/model/ModelFactory.h>
#include <jive/util/utilities.h>
#include <jive/util/FuncUtils.h>
#include <jive/util/Printer.h>

#include "WeakPBCModel2.h"
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
// using jive::geom::BoundaryLine3;
// using jive::geom::Line2;
// using jive::geom::Line3;
// using jive::geom::ShapeFactory;
using jive::model::StateVector;
using jive::util::FuncUtils;


//=========================================================
//    class WeakPBCModel2
//=========================================================
//-----------------------------------------------------------------------
//   static constants
//-----------------------------------------------------------------------

// const char *WeakPBCModel2::COARSEN_METHOD = "coarsenMethod";
const char *WeakPBCModel2::COARSEN_FACTOR = "coarsenFactor";
const char *WeakPBCModel2::NUM_TNODE_PROP = "numTNode";
const char *WeakPBCModel2::ANGLE_PROP = "angle";
const char *WeakPBCModel2::MIN_DIST_PROP = "minDist";
const double WeakPBCModel2::EPS = 1.e-16;
const double WeakPBCModel2::PI = 3.14159265358979323846;

//-----------------------------------------------------------------------
//   constructor
//-----------------------------------------------------------------------

WeakPBCModel2::WeakPBCModel2

	  (const String &name,
	  const Properties &conf,
	  const Properties &props,
	  const Properties &globdat) : PeriodicBCModel::PeriodicBCModel(name, conf, props, globdat)
{
	Properties myProps = props.getProps(myName_);
	Properties myConf = conf.makeProps(myName_);

	myProps.get(cf_, COARSEN_FACTOR);
	myConf.set(COARSEN_FACTOR, cf_);

	myProps.get(numTNode_, NUM_TNODE_PROP);
	myConf.set(NUM_TNODE_PROP, numTNode_);

	angle_ = 0;
	myProps.find(angle_, ANGLE_PROP, -45., 135.);
	myConf.set(ANGLE_PROP, angle_);

	minDist_ = 1.e-6;
	myProps.find(minDist_, MIN_DIST_PROP);
	myConf.set(MIN_DIST_PROP, minDist_);

	nodes_ = XNodeSet::find(globdat);

  // make new NodeSet so as not to mess with nodes_.size() in other models
  // traction dofs will be in the original DofSpace, associated with nodes_
  // but the coordinates of the traction nodes are in tnodes_;

  tnodes_ = jive::fem::newXNodeSet();
  telems_ = jive::fem::newXElementSet(tnodes_);
}

//-----------------------------------------------------------------------
//   Destructor
//-----------------------------------------------------------------------

WeakPBCModel2::~WeakPBCModel2()
{
}

//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------

bool WeakPBCModel2::takeAction(const String &action,
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
  return Super::takeAction(action, params, globdat);
}

//-----------------------------------------------------------------------
//   init_
//-----------------------------------------------------------------------

void WeakPBCModel2::init_(const Properties &globdat)
{
	// Get dofs and cons
	dofs_ = XDofSpace::get(nodes_.getData(), globdat);
	cons_ = Constraints::get(dofs_, globdat);

	dofTypes_.resize(rank_);
  tdofTypes_.resize(rank_);
	box_.resize(rank_ * 2);
	dx0_.resize(rank_);
	dx_.resize(rank_);

	// Add displacement dof types
	dofTypes_[0] = dofs_->addType("dx");
	dofTypes_[1] = dofs_->addType("dy");
	if (rank_ > 2)
		dofTypes_[2] = dofs_->addType("dz");

  // Add traction dof types
  tdofTypes_[0] = dofs_->addType("tx");
  tdofTypes_[1] = dofs_->addType("ty");
  if (rank_ > 2)
    tdofTypes_[2] = dofs_->addType("tz");

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

	// Create boundary element
	String bscheme = "Gauss1";
	bshape_ = BoundaryLine2::getShape("line", bscheme);
	ipCount_ = bshape_->integrationPointCount();
  unodeCount_ = bshape_->nodeCount();
	tnodeCount_ = bshape_->nodeCount();
  udofCount_ = unodeCount_ * rank_;
	tdofCount_ = tnodeCount_ * rank_;
	localrank_ = bshape_->localRank();

	// Create traction mesh
	System::out() << " tnodeIndices_ = [xmin, ymin] = \n";
	findSmallestElement_();
  createTractionMesh_();
  // clearTractionMesh_();
}

//-----------------------------------------------------------------------
//   sortBndNodes_
//-----------------------------------------------------------------------

void WeakPBCModel2::sortBndNodes_()
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
//   findSmallestElement_
//-----------------------------------------------------------------------

void WeakPBCModel2::findSmallestElement_()
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
			dx0_[iy] = (dx < dx0_[iy]) ? dx : dx0_[iy];
		}
	}
}

//-----------------------------------------------------------------------
//   clearTractionMesh_()
//-----------------------------------------------------------------------

void WeakPBCModel2::clearTractionMesh_()
{
  IdxVector idofs(tdofTypes_.size() * tnodes_.size());
  IdxVector inodes(tnodes_.size());

  IdxVector allnodes(iarray(tnodes_.size()));

  idx_t ndof = dofs_->collectDofIndices(idofs, allnodes, tdofTypes_);

  dofs_->eraseDofs(idofs[slice(0, ndof)]);

  telems_.clear();
  tnodes_.clear();
}

//-----------------------------------------------------------------------
//   createTractionMesh_()
//-----------------------------------------------------------------------

void WeakPBCModel2::createTractionMesh_()
{
	Vector coords(rank_); // coordinate vector
  IdxVector connect(tnodeCount_); // connectivity

  // Loop over faces of tnodeIndices_
  for (idx_t ix = 0; ix < rank_; ++ix)
	{
		IdxVector inodes(bndNodes_[2 * ix]);     // bndFace_min
		IdxVector jnodes(bndNodes_[2 * ix + 1]); // bndFace_max
		tnodeIndices_[ix].resize(inodes.size() + jnodes.size());

		// Loop over indices of inodes
		idx_t kn = 0;
		for (idx_t in = 0; in < inodes.size(); ++in)
		{
			nodes_.getNodeCoords(coords, inodes[in]);
			coords[ix] = box_[2 * ix + 1];
      tnodeIndices_[ix][kn++] = tnodes_.addNode(coords);
      // tnodeIndices_[ix][kn++] = nodes_.addNode(coords);
		}

		// Loop over indices of jnodes
		for (idx_t jn = 0; jn < jnodes.size(); ++jn)
		{
			nodes_.getNodeCoords(coords, jnodes[jn]);
      tnodeIndices_[ix][kn++] = tnodes_.addNode(coords);
      // tnodeIndices_[ix][kn++] = nodes_.addNode(coords);
    }

		// Get correct index
		idx_t iy = (ix + 1) % 2;

		// Sorting trFace (works only for 2D)
		sortTrFace_(tnodeIndices_[ix], iy);

		// Coarsen the mesh (works only for 2D)
		coarsenMesh_(tnodeIndices_[ix], iy);

		// Add dofs to trFace
		for (idx_t in = 0; in < tnodeIndices_[ix].size(); ++in)
			dofs_->addDofs(tnodeIndices_[ix][in], tdofTypes_);

		// Print to verify
		System::out() << "tnodeIndices_[" << ix << "] = " << tnodeIndices_[ix] << "\n";
	}
  // IdxVector allnodes(iarray(2 * numTNode_));
  // dofs_->addDofs(allnodes, tdofTypes_);
}

//-----------------------------------------------------------------------
//   createTractionMesh2_()
//-----------------------------------------------------------------------

void WeakPBCModel2::createTractionMesh2_()
{
	Vector coords(rank_);         // coordinate vector
	IdxVector connect(tnodeCount_); // connectivity

	JEM_PRECHECK(nodes_.size() > 2 * numTNode_);

  // Loop over faces of tnodeIndices_
  for (idx_t ix = 0; ix < rank_; ++ix)
	{
		idx_t iy = (ix + 1) % 2;
		idx_t ie0 = ix * (numTNode_ - 1);
		double x0 = box_[2 * ix + 1]; // max x-coord on face ix
		double y0 = box_[2 * iy]; // min y-coord on face ix
		coords[ix] = x0;
		double dy = dx_[iy] / (numTNode_ - 1); // length of traction element
    tnodeIndices_[ix].resize(numTNode_);

    System::out() << " x0 = " << x0 << ", y0 = " << y0 << ", dy = " << dy << "\n";

    // Loop over indices of inodes
    for (idx_t in = 0; in < numTNode_; ++in)
		{
			coords[iy] = y0 + in * dy;
			tnodeIndices_[ix][in] = tnodes_.addNode(coords);
      // tnodeIndices_[ix][in] = nodes_.addNode(coords);
    }

		for (idx_t ie = 0; ie < numTNode_ - 1; ++ie)
		{
			connect[0] = ie0 + ie + ix;
			connect[1] = ie0 + ie + ix + 1;
			telems_.addElement(connect);
		}

    // Add dofs to trFace
    for (idx_t in = 0; in < tnodeIndices_[ix].size(); ++in)
      dofs_->addDofs(tnodeIndices_[ix][in], tdofTypes_);

    // Print to verify
    System::out() << "tnodeIndices_[" << ix << "] = " << tnodeIndices_[ix] << "\n";
  }
	// IdxVector allnodes(iarray(2 * numTNode_));
	// dofs_->addDofs(allnodes, tdofTypes_);
}

//-----------------------------------------------------------------------
//   sortbndFace_
//-----------------------------------------------------------------------

void WeakPBCModel2::sortBndFace_(IdxVector &bndFace, const idx_t &index)
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
      if (c0[index] > c1[index])
      {
        bndFace[jn + 1] = bndFace[jn] + bndFace[jn + 1]; // b = a+b
        bndFace[jn] = bndFace[jn + 1] - bndFace[jn];     // a = (a+b)-a = b
        bndFace[jn + 1] = bndFace[jn + 1] - bndFace[jn]; // b = (a+b)-b = a
      }
    }
  }
}

//-----------------------------------------------------------------------
//   sortTractionFace_
//-----------------------------------------------------------------------

void WeakPBCModel2::sortTrFace_(FlexVector &trFace, const idx_t &index)
{
  Vector c0(rank_); // coordinate vector of node "0"
  Vector c1(rank_); // coordinate vector of node "1"

  // Bubblesort algorithm
  for (idx_t in = 0; in < trFace.size(); ++in)
  {
    for (idx_t jn = 0; jn < trFace.size() - 1 - in; ++jn)
    {
      // Get nodal coordinatess
      tnodes_.getNodeCoords(c1, trFace[jn + 1]);
      tnodes_.getNodeCoords(c0, trFace[jn]);
      // nodes_.getNodeCoords(c1, trFace[jn + 1]);
      // nodes_.getNodeCoords(c0, trFace[jn]);

      // Swap indices if necessary
      if (c0[index] > c1[index])
      {
        trFace[jn + 1] = trFace[jn] + trFace[jn + 1]; // b = a+b
        trFace[jn] = trFace[jn + 1] - trFace[jn];     // a = (a+b)-a = b
        trFace[jn + 1] = trFace[jn + 1] - trFace[jn]; // b = (a+b)-b = a
      }
    }
  }
}

//-----------------------------------------------------------------------
//   coarsenMesh_
//-----------------------------------------------------------------------

void WeakPBCModel2::coarsenMesh_(FlexVector &trFace, const idx_t &index)
{
	using jem::numeric::norm2;

	Vector c0(rank_); // coordinate vector of node "0"
	Vector c1(rank_); // coordinate vector of node "1"
	Vector cn(rank_); // coordinate vector of node "n"
	tnodes_.getNodeCoords(cn, trFace.back());

	// Minimum coarsening factor for Neumann boundary conditions
	double cf_min = (dx0_[0] + dx0_[1]) / rank_ / max(dx_[0], dx_[1]);
	cf_ = (cf_ < cf_min) ? cf_min : cf_;
	// System::out() << "Minimum coarsening cf is " << cf_min << "\n";

	// Tolerance for removing nodes
	double dx = (dx0_[0] + dx0_[1]) / (2 * cf_);

	// Loop over indices of trFace
	int jn = 0;
	for (Iter in = trFace.begin(); in < trFace.end() - 1; ++in)
	{
		// Get nodal coordinates
		tnodes_.getNodeCoords(c0, trFace[jn]);
		tnodes_.getNodeCoords(c1, trFace[jn + 1]);

		// Delete indices until c1 - c0 > dx
		while (norm2(c0 - c1) < min(dx, dx_[index]) && jn < trFace.size())
		{
			// Delete current node index
			trFace.erase(in + 1);
			// Assign next nodal coordinates to c1
			tnodes_.getNodeCoords(c1, trFace[jn + 1]);
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

void WeakPBCModel2::getTractionMeshNodes_(IdxVector &connect,
										 const Vector &x,
										 const idx_t &face)
{
	Matrix coords(rank_, tnodeCount_);

	// Map face onto ix
	int ix = face / 2;

	// Implementation for two dimensions
	if (rank_ == 2)
	{
		// Loop over indices of trFace
		for (idx_t in = 0; in < tnodeIndices_[ix].size() - 1; ++in)
		{
			// get coords of nodes in trFace[in:in+2]
			connect[0] = tnodeIndices_[ix][in];
			connect[1] = tnodeIndices_[ix][in + 1];
			tnodes_.getSomeCoords(coords, connect);

			// Get correct index
			idx_t iy = (ix + 1) % 2;

			// Check if c0[iy] < x[iy] < c1[iy]
			if ((coords(iy, 0) <= x[iy]) && (x[iy] <= coords(iy, 1)))
				return;
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

void WeakPBCModel2::augmentMatrix_(Ref<MatrixBuilder> mbuilder,
								  const Vector &fint,
								  const Vector &solu)
{
	// Variables related to both U mesh and T mesh:
	IdxVector connect(tnodeCount_);    // node indices = [ node1 node2 ]
	Matrix coords(rank_, tnodeCount_); // node coordinates = [ x[1] x[2] ], x = [x y z]'

	// Variables related to element on U mesh:
	IdxVector idofs(udofCount_); // dof indices related to U
	Vector w(ipCount_);         // Integration weights = [ jw[1] jw[2] ]
	Matrix n(unodeCount_, ipCount_);  // shape functions of U mesh: N = [ n[1] n[2] ], n[ip] = [n1 n2]'
	Matrix N(rank_, udofCount_); // N matrix of U mesh [n1 0 n2 0; 0 n1 0 n2]
	Matrix X(rank_, ipCount_);  // global coordinates of IP; X = [ x[1] x[2] ], x[ip] = [x y z]'
	N = 0.0;

	// Variables related to element on T mesh:
	IdxVector jdofs(tdofCount_); // dof indices related to T
	Vector xi(localrank_);  // local coordinates of given X[ip]; u = xi
	Vector h(tnodeCount_);        // shape functions of T; h = [h1 h2]'
	Matrix H(rank_, tdofCount_); // H matrix of T mesh [h1 0 h2 0; 0 h1 0 h2]
	H = 0.0;

	// Matrices and vector to be assembled:
	Vector tr(tdofCount_);         // Vector of tractions of each element
	Vector u(udofCount_);          // Vector of displacements of each element
	Matrix Ke(udofCount_, tdofCount_);  // Ke = w[ip]*N[ip]*H[ip]
	Matrix KeT( Ke.transpose() ); // KeT = w[ip]*H[ip]*N[ip]
	Vector fe(tdofCount_);         // fe = Ke*t or KeT*u

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
			bshape_->getGlobalIntegrationPoints(X, coords);
			dofs_->getDofIndices(idofs, connect, dofTypes_);
			n = bshape_->getShapeFunctions();

			for (idx_t ip = 0; ip < ipCount_; ip++)
			{
				// Get jdofs from T-mesh
				getTractionMeshNodes_(connect, X[ip], face);
				dofs_->getDofIndices(jdofs, connect, tdofTypes_);
				tnodes_.getSomeCoords(coords, connect);

				// Get N-matrix from U-mesh
				N(0, 0) = N(1, 1) = n(0, ip);
				N(0, 2) = N(1, 3) = n(1, ip);

				// Get H-matrix from T-mesh
				int ix = face / 2;
				X(ix, ip) = box_[2 * ix + 1];
				double eps = dx0_[ix] / 100;
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
	Matrix Ht(tdofCount_, rank_);
	Vector u_corner(rank_);
	Vector u_fixed(rank_);
	IdxVector kdofs(rank_);
	dofs_->getDofIndices(kdofs, ifixed_, dofTypes_);
	u_fixed = solu[kdofs];

	// Loop over faces of tnodeIndices_
	for (idx_t ix = 0; ix < rank_; ++ix)
	{
		// Obtain u_corner
		IdxVector idofs(rank_);
		dofs_->getDofIndices(idofs, masters_[ix], dofTypes_);
		u_corner = solu[idofs];

		// Loop over indices of trFace
		for (idx_t jn = 0; jn < tnodeIndices_[ix].size() - 1; ++jn)
		{
			// Get jdofs, w and H from traction mesh
			connect[0] = tnodeIndices_[ix][jn];
			connect[1] = tnodeIndices_[ix][jn + 1];
			tnodes_.getSomeCoords(coords, connect);
			bshape_->getIntegrationWeights(w, coords);
			dofs_->getDofIndices(jdofs, connect, tdofTypes_);
			n = bshape_->getShapeFunctions();

			for (idx_t ip = 0; ip < ipCount_; ip++)
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
//   newWeakPBCModel2
//-----------------------------------------------------------------------

static Ref<Model> newWeakPBCModel2

	(const String &name,
	 const Properties &conf,
	 const Properties &props,
	 const Properties &globdat)

{
	return newInstance<WeakPBCModel2>(name, conf, props, globdat);
}

//-----------------------------------------------------------------------
//   declareWeakPBCModel2
//-----------------------------------------------------------------------

void declareWeakPBCModel2()

{
	using jive::model::ModelFactory;
	ModelFactory::declare("WeakPBC2", &newWeakPBCModel2);
}
