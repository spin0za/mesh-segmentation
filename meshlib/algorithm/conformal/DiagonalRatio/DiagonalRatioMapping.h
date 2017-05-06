/*!
*      \file DiagonalRatioMapping.h
*      \brief Algorithm for Computing Conformal Mappings based on diagonal ratio
*	   \author David Gu
*      \date Document 07/22/2011
*
*		Algorithm for computing Conformal Mappings based on diagonal ratio
*
*/

#ifndef _DIAGONAL_RATIO_MAPPING_H_
#define _DIAGONAL_RATIO_MAPPING_H_

#include  <math.h>
#include <queue>
#include <list>
#include <vector>
#include <complex>
#include "DiagonalRatioMesh.h"
#include "Structure/Structure.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace MeshLib
{

namespace Holomorphy
{

/*! \brief CDiagonalRatioMapping class
 *
 *	Compute Conformal Mappings based on diagonal ratio
 */
template<typename M>
class CDiagonalRatioMapping
{
public:
	/*! CDiagonalRatioMapping constructor
	*	\param mesh, input mesh
	*/
	CDiagonalRatioMapping( M* pMesh );

	/*! CDiagonalRatioMapping destructor
	*/
	~CDiagonalRatioMapping();
	
	/*! Compute the Conformal Mapping based on diagonal ratio
	*/
	void map();

protected:

	/*! Total number of fixed vertices */
	int  m_nfixed_vertices;

	/*! Total number of free vertices */
	int  m_nfree_vertices;

	/*! Total number of interior edges */
	int m_ninterior_edges;

protected:
	/*!	Pointer to the input mesh
	 */
	M * m_pMesh;

	/*!	Compute the diagonal ratio for one face
	*/
	std::complex<double> _diagonal_ratio( typename M::CEdge * pE );
	/*!	write the uv to the vertex string
	 */
	void _write_uv( M * pMesh );
	/*!	label vertex indices
	 */
	void _label_vertex_edge_index();

	/*!	set fixed vertices' uv coordinates
	 */
	void _set_fixed_vertices();
};

//CDiagonalRatioMapping constructor
//\param mesh is the input mesh
template<typename M>
CDiagonalRatioMapping<M>::CDiagonalRatioMapping( M* pMesh )
{
	m_pMesh = pMesh;


	COperator<M> stru( m_pMesh );
	stru._embedding_2_metric();
	stru._metric_2_diagonal_ratio();
};

//CDiagonalRatioMapping destructor
template<typename M>
CDiagonalRatioMapping<M>::~CDiagonalRatioMapping()
{
};


// Compute holomorphic 1-forms
template<typename M>
void CDiagonalRatioMapping<M>::map()
{
	//label vertex and edge index
	_label_vertex_edge_index();

	//set fixed vertices
	_set_fixed_vertices();

	std::vector<Eigen::Triplet<std::complex<double> > >  M_coefficients;
	std::vector<Eigen::Triplet<std::complex<double> > > MT_coefficients;
	std::vector<Eigen::Triplet<std::complex<double> > >  F_coefficients;
	
	//construct matrix
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * pE = *eiter;

		//only consider interior edges
		
		if( m_pMesh->isBoundary( pE ) ) continue;

		M::CHalfEdge * pH = m_pMesh->edgeHalfedge( pE, 0 );
		M::CHalfEdge * pN = m_pMesh->faceNextCcwHalfEdge( pH );
		M::CHalfEdge * pD = m_pMesh->edgeHalfedge( pE, 1 );
		M::CHalfEdge * pW = m_pMesh->faceNextCcwHalfEdge( pD );
		
		M::CVertex * pS = m_pMesh->halfedgeSource( pH );
		M::CVertex * pT = m_pMesh->halfedgeTarget( pH );
		M::CVertex * pL = m_pMesh->halfedgeTarget( pN );
		M::CVertex * pR = m_pMesh->halfedgeTarget( pW );

	
		std::complex<double> rho = pE->rho();

		if( !pT->fixed() )
		{
			 M_coefficients.push_back( Eigen::Triplet<std::complex<double> >( pE->idx(), pT->idx(), 1.0) );
			MT_coefficients.push_back( Eigen::Triplet<std::complex<double> >( pT->idx(), pE->idx(), 1.0) );
		}
			//M.AddElementTail( pE->idx(), pT->idx(), +1 );
		else
		{
			F_coefficients.push_back( Eigen::Triplet<std::complex<double> >( pE->idx(), pT->idx(), 1.0) );
		}
			//F.AddElementTail( pE->idx(), pT->idx(), +1 );

		if( !pS->fixed() )
		{
			 M_coefficients.push_back( Eigen::Triplet<std::complex<double> >( pE->idx(), pS->idx(), -1.0) );
			MT_coefficients.push_back( Eigen::Triplet<std::complex<double> >( pS->idx(), pE->idx(), -1.0) );
		}
			//M.AddElementTail( pE->idx(), pS->idx(), -1 );
		else
		{
			F_coefficients.push_back( Eigen::Triplet<std::complex<double> >( pE->idx(), pS->idx(), -1.0) );
		}
			//F.AddElementTail( pE->idx(), pS->idx(), -1 );

		if( !pL->fixed() )
		{
			 M_coefficients.push_back( Eigen::Triplet<std::complex<double> >( pE->idx(), pL->idx(), -rho) );
			 MT_coefficients.push_back( Eigen::Triplet<std::complex<double> >( pL->idx(), pE->idx(), std::conj(-rho)) );
		}

			//M.AddElementTail( pE->idx(), pL->idx(), -rho );
		else
		{
			F_coefficients.push_back( Eigen::Triplet<std::complex<double> >( pE->idx(), pL->idx(), -rho) );
		}
			//F.AddElementTail( pE->idx(), pL->idx(), -rho );

		if( !pR->fixed() )
		{
			M_coefficients.push_back( Eigen::Triplet<std::complex<double> >( pE->idx(), pR->idx(), rho) );
		   MT_coefficients.push_back( Eigen::Triplet<std::complex<double> >( pR->idx(), pE->idx(), std::conj(rho)) );
		}

			//M.AddElementTail( pE->idx(), pR->idx(),  rho );
		else
		{
			F_coefficients.push_back( Eigen::Triplet<std::complex<double> >( pE->idx(), pR->idx(), rho) );
		}
			//F.AddElementTail( pE->idx(), pR->idx(),  rho );
	};

	//construct sparse matrices
	Eigen::SparseMatrix<std::complex<double>> M( m_ninterior_edges, m_nfree_vertices );
	M.setZero();
	//construct sparse matrices
	Eigen::SparseMatrix<std::complex<double>> MT( m_nfree_vertices, m_ninterior_edges );
	MT.setZero();

	Eigen::SparseMatrix<std::complex<double>> F( m_ninterior_edges, m_nfixed_vertices );
	F.setZero();
	
	 M.setFromTriplets( M_coefficients.begin(), M_coefficients.end());
	MT.setFromTriplets(MT_coefficients.begin(), MT_coefficients.end());
	 F.setFromTriplets(F_coefficients.begin(), F_coefficients.end());

	Eigen::VectorXcd b( m_nfixed_vertices );

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		if( !pV->fixed() ) continue;
		b(pV->idx() ) = std::complex<double>(-pV->huv()[0], -pV->huv()[1]);
	}

	


	Eigen::SparseMatrix<std::complex<double> > P = MT*M;
	Eigen::VectorXcd B = MT*F*b;	

	std::cerr << "Start factoring" << std::endl;
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<std::complex<double>> > solver;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>> > solver;
	solver.compute(P);
	std::cerr << "End factoring" << std::endl;

	if(solver.info()!=Eigen::Success ) {
		std::cerr << "Error in factoring" << std::endl;
		return;
	}

	std::cerr << "Start Solving" << std::endl;
	Eigen::VectorXcd x = solver.solve(B);
	std::cerr << "End Solving" << std::endl;

	if(solver.info()!=Eigen::Success ) {
		std::cerr << "Error in factoring" << std::endl;
		return;
	}

/*
	for( size_t i = 0; i < x.size(); i ++ )
	{
		std::cout << x[i] << " ";
	}
*/
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		if( pV->fixed() ) continue;
		std::complex<double> c = x(pV->idx());
		pV->huv() = CPoint2( c.real(), c.imag() );
	}

	_write_uv( m_pMesh );

};


template<typename M>
void CDiagonalRatioMapping<M>::_label_vertex_edge_index()
{
	//find all fixed vertices
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;

		int valence = 0;

		for( M::VertexEdgeIterator veiter( pV ); !veiter.end(); veiter ++ )
		{
			M::CEdge * pE = *veiter;
			if( pE->sharp() )
				valence ++;
		}

		pV->fixed() = (valence >= 2 );
	}

	//label vertex index
	int fixed_vid = 0;
	int free_vid  = 0;

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		if( pV->fixed() ) 
		{
			pV->idx() = fixed_vid ++;
		}
		else
		{
			pV->idx() = free_vid  ++;
		}
	}
	//total number of fixed vertices
	m_nfixed_vertices = fixed_vid;

	//total number of free vertices
	m_nfree_vertices  = free_vid;

	//label edge index
	int eid = 0;	
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * pE = *eiter;
		if( m_pMesh->isBoundary( pE ) ) continue;
		pE->idx() = eid ++;
	}

	//total number of edges
	m_ninterior_edges = eid;
};


template<typename M>
void CDiagonalRatioMapping<M>::_write_uv( M * pMesh )
{

	for( M::MeshVertexIterator viter( pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		std::string &m_string = pV->string();

		CParser parser( m_string );
		parser._removeToken( "uv" );

		parser._toString( m_string );
		
		std::string line;
		std::stringstream iss(line);
		iss << "uv=(" << pV->huv()[0] << " " << pV->huv()[1] << ")";

		if( m_string.length() > 0 )
		{
			m_string += " ";
			m_string += iss.str();
		}
		else
		{
			m_string = iss.str();
		}
	}
};

template<typename M>
void CDiagonalRatioMapping<M>::_set_fixed_vertices()
{
	std::vector<M::CVertex*> fixed_vertices;

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		if( pV->fixed()) fixed_vertices.push_back( pV );
	}

	fixed_vertices[0]->huv() = CPoint2(0, 0 );
	fixed_vertices[1]->huv() = CPoint2(0, 1 );

};

} //namespace Holomorphy
} //namespace MeshLib
#endif