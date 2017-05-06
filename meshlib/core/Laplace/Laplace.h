/*!
*      \file Laplace.h
*      \brief Solve Laplace equation with Dirichlet boundary condition
*	   \author David Gu
*      \date Document 12/08/2013
*
*/

#ifndef _LAPLACE_H_
#define _LAPLACE_H_

#include <map>
#include <vector>
#include <complex>

#include "Mesh/BaseMesh.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "mesh/boundary.h"
#include "Parser/parser.h"
#include <Eigen/Sparse>


namespace MeshLib
{
	template<typename M>
	class CLaplace
	{
		public:

		/*! \brief COperator constructor
		 *  \param pMesh the input mesh
		 */
		  CLaplace( M * pMesh ) { m_pMesh = pMesh; };
		/*! \brief COperator destructor
		 */
		  ~CLaplace(){};
		
		  /*!	\brief Solve Laplace Equation
		   *    \param input : the edge weight is set, the fixed vertex u is set
		   *    \param output: the vertex u is set
		   */
		  void solve();

		protected:
			M * m_pMesh;
	};

/*!	CLaplace constructor 
*	Count the number of interior vertices, boundary vertices and the edge weight
*
*/
template<typename M>
void CLaplace<M>::solve()
{
	int vid  = 0; //interior vertex ID 
	int vfid = 0; //boundary vertex ID

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * pV = *viter;
		if( pV->fixed() )
		{
			pV->idx() = vfid ++;
		}
		else
		{
			pV->idx() = vid  ++;
		}
	}

	int free_vertices = vid;
	int fixed_vertices = vfid;
	


	std::vector<Eigen::Triplet<double> > A_coefficients;
	std::vector<Eigen::Triplet<double> > B_coefficients;

	
	//set the matrix A
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); ++ eiter )
	{
		M::CEdge * pE = * eiter;

		M::CVertex * pV = m_pMesh->edgeVertex1( pE );
		M::CVertex * pW = m_pMesh->edgeVertex2( pE );

		if( pV->fixed() && pW->fixed() ) continue;
		
		if( pV->fixed() && !pW->fixed() )
		{
			M::CVertex * pT = pV;
			pV = pW;
			pW = pT;
		}
		
		int vid = pV->idx();
		int wid = pW->idx();

		double w = pE->weight();

		if( pW->fixed() )
		{
			Eigen::Triplet<double> e(vid,wid,w);
			B_coefficients.push_back( e );
		}
		else
		{
			A_coefficients.push_back( Eigen::Triplet<double>(vid,wid, -w) );
			A_coefficients.push_back( Eigen::Triplet<double>(wid,vid, -w) );
		}
	};

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * pV = *viter;
		if( pV->fixed() ) continue;
		int vid = pV->idx();

		double sw = 0;
		for( M::VertexEdgeIterator eiter( pV ); !eiter.end(); ++ eiter )
		{
			M::CEdge * pE = *eiter;
			double w = pE->weight();
			sw += w;
		}
		A_coefficients.push_back( Eigen::Triplet<double>(vid,vid, sw ) );
	}


	Eigen::SparseMatrix<double> A( free_vertices, free_vertices );
	A.setZero();

	Eigen::SparseMatrix<double> B( free_vertices, fixed_vertices );
	B.setZero();
	A.setFromTriplets(A_coefficients.begin(), A_coefficients.end());
	B.setFromTriplets(B_coefficients.begin(), B_coefficients.end());


	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	std::cerr << "Eigen Decomposition" << std::endl;
	solver.compute(A);
	std::cerr << "Eigen Decomposition Finished" << std::endl;
	
	if( solver.info() != Eigen::Success )
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}


	Eigen::VectorXd b(fixed_vertices);

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * pV = *viter;
		if( !pV->fixed() ) continue;
		int id = pV->idx();
		b(id) = pV->u();
	}

	Eigen::VectorXd c(free_vertices);
	c = B * b;

	Eigen::VectorXd x = solver.solve(c);
	if( solver.info() != Eigen::Success )
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}

	//set the images of the harmonic map to interior vertices
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * pV = *viter;
		if( pV->fixed() ) continue;
		int id = pV->idx();
		pV->u() = x(id);
	}
};


}; //end meshLib

#endif _LAPLACE_H_