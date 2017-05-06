/*!
*      \file DiagonalRatioSymmetryMapping.h
*      \brief Algorithm for Computing Symmetric Conformal Mappings based on diagonal ratio
*	   \author David Gu
*      \date Document 11/01/2013
*
*		Algorithm for computing Symmetric Conformal Mappings based on diagonal ratio
*
*/

#ifndef _DIAGONAL_RATIO_SYMMETRY_MAPPING_H_
#define _DIAGONAL_RATIO_SYMMETRY_MAPPING_H_

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

/*! \brief CDiagonalRatioSymmetryMapping class
 *
 *	Compute Symmetric Conformal Mappings based on diagonal ratio
 */
template<typename M>
class CDiagonalRatioSymmetryMapping
{
public:
	/*! CDiagonalRatioSymmetryMapping constructor
	*	\param mesh, input mesh
	*/
	CDiagonalRatioSymmetryMapping( M* pMesh ){ m_pMesh = pMesh; };

	/*! CDiagonalRatioMapping destructor
	*/
	~CDiagonalRatioSymmetryMapping(){};
	
	/*! Compute the Conformal Mapping based on diagonal ratio
	*/
	void map();
	/*! input the marker file
	 */
	void _input_markers( const char * marker_file );

protected:

	/*! Total number of free vertices */
	int  m_nvertices;

	/*! Total number of interior edges */
	int m_ninterior_edges;

	/*! number of vertices on the axis */
	int m_nVerts_on_axis;
	std::vector<typename M::CVertex*> m_verts_on_axis;

	/*! number of symmetric vertex pairs */
	int m_nVerts_pairs;
	std::vector<std::pair<typename M::CVertex*, typename M::CVertex*> > m_vert_pairs;

	/*! number of anchors */
	int m_nVerts_anchors;
	std::vector<std::pair<typename M::CVertex*, CPoint2> > m_vert_anchors;
	
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

	/*! build the linear system
	 */
	void _form_linear_system();

};

template<typename M>
void CDiagonalRatioSymmetryMapping<M>::_write_uv( M * pMesh )
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
void CDiagonalRatioSymmetryMapping<M>::_label_vertex_edge_index()
{
	int vid = 0;
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		pV->idx() = vid  ++;
	}
	m_nvertices = vid;

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
void CDiagonalRatioSymmetryMapping<M>::_input_markers( const char * marker_file )
{
	std::fstream is( marker_file, std::fstream::in );

	if( is.fail() )
	{
		std::cerr << "Error in opening file " << marker_file << std::endl;
		return;
	}



	is >> m_nVerts_on_axis >> m_nVerts_pairs >> m_nVerts_anchors;
	for( int i = 0; i < m_nVerts_on_axis; i ++ )
	{
		int vid;
		is >> vid;
		M::CVertex * pV = m_pMesh->idVertex( vid );
		m_verts_on_axis.push_back( pV );
	}

	for( int i = 0; i < m_nVerts_pairs; i ++ )
	{
		int vid;
		is >> vid;
		M::CVertex * pV = m_pMesh->idVertex( vid );
		int wid;
		is >> wid;
		M::CVertex * pW = m_pMesh->idVertex( wid );

		m_vert_pairs.push_back( std::pair<M::CVertex*,M::CVertex*>( pV, pW ) );
	}

	for( int i = 0; i < m_nVerts_anchors; i ++ )
	{
		int vid;
		is >> vid;
		M::CVertex * pV = m_pMesh->idVertex( vid );
		
		double x, y;
		is >> x >> y;
		m_vert_anchors.push_back( std::pair<M::CVertex*,CPoint2>( pV, CPoint2(x,y) ) );
	}
};


template<typename M>
void CDiagonalRatioSymmetryMapping<M>::_form_linear_system()
{
	std::vector<Eigen::Triplet<std::complex<double> > >  M_coefficients;

	//construct matrix
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * pE = *eiter;

		//only consider interior edges
		
		if( m_pMesh->isBoundary( pE ) ) continue;

		CHalfEdge * pH = m_pMesh->edgeHalfedge( pE, 0 );
		CHalfEdge * pN = m_pMesh->faceNextCcwHalfEdge( pH );
		CHalfEdge * pD = m_pMesh->edgeHalfedge( pE, 1 );
		CHalfEdge * pW = m_pMesh->faceNextCcwHalfEdge( pD );
		
		CDRVertex * pS = m_pMesh->halfedgeSource( pH );
		CDRVertex * pT = m_pMesh->halfedgeTarget( pH );
		CDRVertex * pL = m_pMesh->halfedgeTarget( pN );
		CDRVertex * pR = m_pMesh->halfedgeTarget( pW );

	
		std::complex<double> rho = pE->rho();

		M_coefficients.push_back( Eigen::Triplet<std::complex<double> >( pE->idx(), pT->idx(), 1.0)  );
		M_coefficients.push_back( Eigen::Triplet<std::complex<double> >( pE->idx(), pS->idx(), -1.0) );
		M_coefficients.push_back( Eigen::Triplet<std::complex<double> >( pE->idx(), pL->idx(), -rho) );
		M_coefficients.push_back( Eigen::Triplet<std::complex<double> >( pE->idx(), pR->idx(), rho)  );
	};

	std::vector<Eigen::Triplet<double> > entries;

	for( size_t i = 0; i < M_coefficients.size(); i ++ )
	{
		Eigen::Triplet<std::complex<double> > & e = M_coefficients[i];
		entries.push_back( Eigen::Triplet<double>( e.row(), e.col(), e.value().real() ) );
		entries.push_back( Eigen::Triplet<double>( e.row() + m_ninterior_edges, e.col()+ m_nvertices, e.value().real() ) );
		entries.push_back( Eigen::Triplet<double>( e.row() + m_ninterior_edges, e.col(),  e.value().imag() ) );
		entries.push_back( Eigen::Triplet<double>( e.row(), e.col() + m_nvertices, -e.value().imag() ) );
	}

	std::vector<Eigen::Triplet<double>> & axis = entries;
	//symmetry axis
	int row_id = 2 * m_ninterior_edges;
	for( size_t i = 0; i < m_verts_on_axis.size(); i ++ )
	{
		M::CVertex * pV = m_verts_on_axis[i];
		axis.push_back( Eigen::Triplet<double>( row_id ++, pV->idx(), 1.0 ) );
	}

	//symmetry markers
	std::vector<Eigen::Triplet<double>> & dual = entries;

	for( size_t i = 0; i < m_vert_pairs.size(); i ++ )
	{
		std::pair<M::CVertex*,M::CVertex*> & p = m_vert_pairs[i];

		M::CVertex * pV1 = p.first;
		M::CVertex * pV2 = p.second;

		dual.push_back( Eigen::Triplet<double>( row_id   , pV1->idx(),  1.0 ) );
		dual.push_back( Eigen::Triplet<double>( row_id ++, pV2->idx(),  1.0 ) );

		dual.push_back( Eigen::Triplet<double>( row_id   , pV1->idx() + m_nvertices,  1.0 ) );
		dual.push_back( Eigen::Triplet<double>( row_id ++, pV2->idx() + m_nvertices, -1.0 ) );
	}
	//markers
	std::vector<Eigen::Triplet<double>> & anchors = entries;

	for( size_t i = 0; i < m_vert_anchors.size(); i ++ )
	{
		std::pair<M::CVertex*, CPoint2> & p = m_vert_anchors[i];
		M::CVertex * pV = p.first;

		anchors.push_back( Eigen::Triplet<double>( row_id ++, pV->idx(), 1.0 ) ); //real part
		anchors.push_back( Eigen::Triplet<double>( row_id ++, pV->idx() + m_nvertices, 1.0 ) ); //imaginary part
	}
/*
	std::vector<Eigen::Triplet<double> > total_entries;
	total_entries.reserve( entries.size() + axis.size() + dual.size() + anchors.size() );

	total_entries.insert( entries.begin(), entries.end() );
	total_entries.insert( axis.begin(), axis.end() );
	total_entries.insert( dual.begin(), dual.end() );
	total_entries.insert( anchors.begin(), anchors.end() );
*/
	int rows = m_ninterior_edges * 2 + m_nVerts_on_axis + m_nVerts_pairs * 2 + m_nVerts_anchors *2;
	int cols = m_nvertices * 2;

	Eigen::SparseMatrix<double> M( rows, cols );
	M.setZero();
	//construct sparse matrices
	Eigen::SparseMatrix<double> MT( cols, rows );
	MT.setZero();

	//M.setFromTriplets( total_entries.begin(), total_entries.end());
	M.setFromTriplets( entries.begin(), entries.end());
	MT = Eigen::SparseMatrix<double>( M.transpose() );


	Eigen::VectorXd b( rows );
	b.setZero();

	int n = rows - m_nVerts_anchors * 2;

	for( size_t i = 0; i < m_vert_anchors.size(); i ++ )
	{
		std::pair<M::CVertex*,CPoint2> & p = m_vert_anchors[i];
		CPoint2 pt = p.second;
		b( n + i * 2 + 0 ) = pt[0];
		b( n + i * 2 + 1 ) = pt[1];
	}

	Eigen::VectorXd B( cols );
	B.setZero();
	B = MT * b;

	Eigen::SparseMatrix<double> P = MT*M;

	std::cerr << "Start factoring" << std::endl;
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
	solver.compute(P);
	std::cerr << "End factoring" << std::endl;

	if(solver.info()!=Eigen::Success ) {
		std::cerr << "Error in factoring" << std::endl;
		return;
	}

	std::cerr << "Start Solving" << std::endl;
	Eigen::VectorXd x = solver.solve(B);
	std::cerr << "End Solving" << std::endl;

	if(solver.info()!=Eigen::Success ) {
		std::cerr << "Error in factoring" << std::endl;
		return;
	}

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		int d = pV->idx();
		pV->huv() = CPoint2( x(d), x(d+m_nvertices) );
	}

};


/*! Assume the markers have been input already
 */
template<typename M>
void CDiagonalRatioSymmetryMapping<M>::map()
{

	CStructure<M, M::CVertex, M::CEdge, M::CFace, M::CHalfEdge> stru( m_pMesh );
	stru._embedding_2_metric();
	stru._metric_2_diagonal_ratio();
	//label vertex and edge index
	_label_vertex_edge_index();

	_form_linear_system();
	_write_uv( m_pMesh );

};




} //namespace Holomorphy
} //namespace MeshLib

#endif