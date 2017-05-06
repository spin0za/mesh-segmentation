/*!
*      \file DiagonalRatioHolomorphicForm.h
*      \brief Algorithm for computing holomorphic 1-forms based on diagonal ratio
*	   \author David Gu
*      \date Document 06/29/2011
*
*		Algorithm for computing holomorphic 1-forms based on diagonal ratio
*
*/

#ifndef _DIAGONAL_RATIO_HOLOMORPHIC_FORM_H_
#define _DIAGONAL_RATIO_HOLOMORPHIC_FORM_H_

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

/*! \brief CDiagonalRatioHolomorphicForm class
 *
 *	Compute holomorphic forms based on diagonal ratio
 */
template<typename M>
class CDiagonalRatioHolomorphicForm
{
public:
	/*! CDiagonalRatioHolomorphicForm constructor
	*	\param meshes the list of meshes with stores the harmonic 1-forms, which form the basis 
	*   of the first cohomology group of the mesh
	*/
	CDiagonalRatioHolomorphicForm( std::list<M*> & meshes );
	/*! CDiagonalRatioHolomorphicForm destructor
	*/
	~CDiagonalRatioHolomorphicForm();
	/*!	The list of meshes storing harmonic form bases
	*/
	std::vector<M*> & meshes() { return m_meshes; };

	/*! Compute the holomorphic 1-form based on diagonal ratio
	*/
	void conjugate();

protected:
	/*!	The list of meshes storing harmonic form bases
	*/
	std::vector<M*> m_meshes;
	/*!	Pointer to the current mesh
	 */
	M * m_pMesh;

	/*!	Compute the diagonal ratio for one face
	*/
	std::complex<double> _diagonal_ratio( typename M::CEdge * pE );
	/*!	write the duv to the vertex string
	 */
	void _write_duv( M * pMesh );

private:

/*! 
 *	integrate 1-form on a chain 
 */
	double _integrate( M * form, std::vector<typename M::CVertex*> chain );

/*!	Evalulate a 1-form on a halfedge
 */
double _evaluate( M * pMesh, int source_vert_id, int target_vert_id );
};


//CDiagonalRatioHolomorphicForm constructor
//\param meshes are the basis of harmonic 1-form group
template<typename M>
CDiagonalRatioHolomorphicForm<M>::CDiagonalRatioHolomorphicForm( std::list<M*> & meshes )
{

	for( std::list<M*>::iterator miter = meshes.begin(); miter != meshes.end(); miter ++ )
	{
		M * pM = *miter;
		m_meshes.push_back( pM );
	}

	m_pMesh = m_meshes[0];

	//label vertex index
	int vid = 0;	
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		pV->idx() = vid ++;
	}

	//label edge index
	int eid = 0;	
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * pE = *eiter;
		pE->idx() = eid ++;
	}

	COperator<M> stru( m_pMesh );
	stru._embedding_2_metric();
	stru._metric_2_diagonal_ratio();
};

//CDiagonalRatioHolomorphicForm destructor
template<typename M>
CDiagonalRatioHolomorphicForm<M>::~CDiagonalRatioHolomorphicForm()
{
};


// Compute holomorphic 1-forms
template<typename M>
void CDiagonalRatioHolomorphicForm<M>::conjugate()
{
	size_t nM = m_meshes.size();

	for( size_t k = 0; k < m_meshes.size(); k ++ )
	{
		M * pM = m_meshes[k];

		std::vector<Eigen::Triplet<std::complex<double>>>  M_coefficients;
		std::vector<Eigen::Triplet<std::complex<double>>> MT_coefficients;

		Eigen::VectorXcd b( pM->numEdges() );
		int id = 0;
		
		
		//construct matrix
		for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
		{
			M::CEdge * pE = *eiter;			

			M::CHalfEdge * pH = m_pMesh->edgeHalfedge( pE, 0 );
			M::CHalfEdge * pN = m_pMesh->faceNextCcwHalfEdge( pH );
			M::CHalfEdge * pD = m_pMesh->edgeHalfedge( pE, 1 );
			M::CHalfEdge * pW = m_pMesh->faceNextCcwHalfEdge( pD );
			
			M::CVertex * pS = m_pMesh->halfedgeSource( pH );
			M::CVertex * pT = m_pMesh->halfedgeTarget( pH );
			M::CVertex * pL = m_pMesh->halfedgeTarget( pN );
			M::CVertex * pR = m_pMesh->halfedgeTarget( pW );

		
			std::complex<double> rho = pE->rho();

			M_coefficients.push_back( Eigen::Triplet<std::complex<double>> ( pE->idx(), pT->idx(), +1 ) );
			M_coefficients.push_back( Eigen::Triplet<std::complex<double>> ( pE->idx(), pS->idx(), -1 ) );
			M_coefficients.push_back( Eigen::Triplet<std::complex<double>> ( pE->idx(), pL->idx(), -rho ) );
			M_coefficients.push_back( Eigen::Triplet<std::complex<double>> ( pE->idx(), pR->idx(), +rho ) );

			MT_coefficients.push_back( Eigen::Triplet<std::complex<double>> ( pT->idx(), pE->idx(), +1 ) );
			MT_coefficients.push_back( Eigen::Triplet<std::complex<double>> ( pS->idx(), pE->idx(), -1 ) );
			MT_coefficients.push_back( Eigen::Triplet<std::complex<double>> ( pL->idx(), pE->idx(), std::conj(-rho) ) );
			MT_coefficients.push_back( Eigen::Triplet<std::complex<double>> ( pR->idx(), pE->idx(), std::conj(+rho) ) );
		};



		//construct matrix
		for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
		{
			M::CEdge * pE = *eiter;			

			M::CHalfEdge * pH = m_pMesh->edgeHalfedge( pE, 0 );
			M::CHalfEdge * pN = m_pMesh->faceNextCcwHalfEdge( pH );
			M::CHalfEdge * pD = m_pMesh->edgeHalfedge( pE, 1 );
			M::CHalfEdge * pW = m_pMesh->faceNextCcwHalfEdge( pD );
			
			M::CVertex * pS = m_pMesh->halfedgeSource( pH );
			M::CVertex * pT = m_pMesh->halfedgeTarget( pH );
			M::CVertex * pL = m_pMesh->halfedgeTarget( pN );
			M::CVertex * pR = m_pMesh->halfedgeTarget( pW );
			
			std::complex<double> rho = pE->rho();
		
			std::vector<M::CVertex*> d0;
			std::vector<M::CVertex*> d1;
			
			d0.push_back( pS );
			d0.push_back( pT );

			d1.push_back( pR );
			d1.push_back( pT );
			d1.push_back( pL );

			std::vector< std::complex<double> > lambda;

			for( size_t j = 0; j < m_meshes.size(); j ++ )
			{
				std::complex<double> t = _integrate( m_meshes[j], d0 ) - rho * _integrate( m_meshes[j], d1 );
				lambda.push_back( t );
			}
			

			//b.push_back( -lambda[k] );
			b( id ++ ) = -lambda[k];

			std::vector< std::complex<double> > coefficients;

			for( size_t j = 0; j < m_meshes.size(); j ++ )
			{
				if( j == k ) continue;
				coefficients.push_back( lambda[j] );
			}
			for( size_t j = 0; j < coefficients.size(); j ++ )
			{
				 M_coefficients.push_back( Eigen::Triplet<std::complex<double>> ( pE->idx(), m_pMesh->numVertices()+j, coefficients[j] ) );
				 MT_coefficients.push_back( Eigen::Triplet<std::complex<double>> ( m_pMesh->numVertices()+j, pE->idx(), std::conj(coefficients[j]) ) );
				//M.AddElementTail( pE->idx(), m_pMesh->numVertices()+j, coefficients[j] );
			}
		}


		//construct sparse matrices
		Eigen::SparseMatrix<std::complex<double>> M( pM->numEdges(),  pM->numVertices() + nM - 1 );
		M.setZero();
		//construct sparse matrices
		Eigen::SparseMatrix<std::complex<double>> MT( pM->numVertices() + nM - 1, pM->numEdges() );
		MT.setZero();

		M.setFromTriplets( M_coefficients.begin(), M_coefficients.end());
		MT.setFromTriplets(MT_coefficients.begin(), MT_coefficients.end());
		
		Eigen::SparseMatrix<std::complex<double>> P = MT * M;

		std::cerr << "Start factoring" << std::endl;
		//Eigen::ConjugateGradient<Eigen::SparseMatrix<std::complex<double>> > solver;
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>> > solver;
		solver.compute(P);
		std::cerr << "End factoring" << std::endl;

		if(solver.info()!=Eigen::Success ) {
			std::cerr << "Error in factoring" << std::endl;
			return;
		}
		
		Eigen::VectorXcd B = MT*b;

		std::cerr << "Start Solving" << std::endl;
		Eigen::VectorXcd x = solver.solve(B);
		std::cerr << "End Solving" << std::endl;

		//for( size_t i = 0; i < x.size(); i ++ )
		//{
		//	std::cout << x(i) << " ";
		//}

		std::vector< std::complex<double> > coefficients;
		for( size_t j = 0; j < m_meshes.size() - 1 ; j ++ )
		{
			coefficients.push_back( x(m_pMesh->numVertices()+j) );
		}

		std::vector< std::complex<double> > lambda;

		for( size_t j = 0; j < k ; j ++ )
		{
			lambda.push_back( coefficients[j] );
		}
		lambda.push_back( std::complex<double>(1,0));
		
		for( size_t j = k ; j < coefficients.size(); j ++ )
		{
			lambda.push_back( coefficients[j] );
		}

		for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
		{
			M::CEdge * pE = *eiter;

			M::CVertex * v1 = pM->edgeVertex1( pE );
			M::CVertex * v2 = pM->edgeVertex2( pE );

			std::complex<double> duv(0,0);

			for( size_t i = 0; i < m_meshes.size(); i ++ )
			{
				
				M::CVertex * pW1 = m_meshes[i]->idVertex( v1->id() );
				M::CVertex * pW2 = m_meshes[i]->idVertex( v2->id() );

				M::CEdge * pWE = m_meshes[i]->vertexEdge( pW1, pW2);

				double d = ( m_meshes[i]->edgeVertex1( pWE ) == pW1 )? pWE->du() : -pWE->du();
				duv += lambda[i] * d;
			}
			
			M::CVertex * w1 = pM->idVertex( v1->id() );
			M::CVertex * w2 = pM->idVertex( v2->id() );

			M::CEdge * pWE = pM->vertexEdge( w1, w2 );

			pWE->duv() = duv;
			pWE->duv() += x(v2->idx()) - x(v1->idx());

		}

		_write_duv( pM );

		//break;
	}
};

template<typename M>
void CDiagonalRatioHolomorphicForm<M>::_write_duv( M * pMesh )
{

	for( M::MeshEdgeIterator eiter( pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * pE = *eiter;
		std::string &m_string = pE->string();

		CParser parser( m_string );
		parser._removeToken( "duv" );

		parser._toString( m_string );
		
		std::string line;
		std::stringstream iss(line);
		iss << "duv=(" << pE->duv().real() << " " << pE->duv().imag() << ")";

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
}

template<typename M>
double CDiagonalRatioHolomorphicForm<M>::_integrate( M * form, std::vector<typename M::CVertex*> chain )
{
	double sum = 0; 
	for( size_t i = 0; i < chain.size() - 1 ; i ++ )
	{
		sum += _evaluate( form, chain[i]->id(), chain[i+1]->id() );
	}
	return sum;
}

template<typename M>
double CDiagonalRatioHolomorphicForm<M>::_evaluate(M *pMesh, int source_vert_id, int target_vert_id )
{
	M::CVertex * pS = pMesh->idVertex( source_vert_id );
	M::CVertex * pT = pMesh->idVertex( target_vert_id );
	assert( pS );
	assert( pT );
	M::CEdge * pE = pMesh->vertexEdge( pS, pT );
	assert( pE );
	double du = pE->du();

	if( pMesh->edgeVertex1( pE ) == pS ) return du;
	return -du;
};






} //namespace Holomorphy
} //namespace MeshLib
#endif