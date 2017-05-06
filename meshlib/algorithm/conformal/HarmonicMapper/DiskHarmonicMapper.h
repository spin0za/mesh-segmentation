/*!
*      \file HarmonicMapper.h
*      \brief Algorithm for harmonic mapping
*	   \author David Gu
*      \date Document 10/07/2010
*
*		Simple harmonic map, that maps a topological disk to the unit disk 
*		with minimal harmonic energy.
*/

// Simple harmonic map, that maps topological disk to a unit disk using harmonic map.

#ifndef _DISK_HARMONIC_MAPPER_H_
#define _DISK_HARMONIC_MAPPER_H_

#include <vector>
#include "HarmonicMapperMesh.h"
#include "Operator/Operator.h"

namespace MeshLib
{

namespace Holomorphy
{

/*!
 *	\brief CDiskHarmonicMapper class
 *
 *  Compute the harmonic map by solving Dirichlet problem
 * \f[
 *		\left\{
 *		\begin{array}{ccc}
 *		 \Delta u &\equiv& 0\\
 *		 u|_{\partial \Omega} &=& f
 *		 \end{array}
 *		\right.
 * \f]
 */
	template<typename M>
	class CDiskHarmonicMapper
	{
	public:
		/*!	CHarmonicMapper constructor
		 *	\param pMesh the input mesh
		 */
		CDiskHarmonicMapper(M* pMesh);
		/*!	CHarmonicMapper destructor
		 */
		~CDiskHarmonicMapper();
		/*!  Compute the harmonic map using direct method
		 */
		void _map();
		/*!	Iterative method compute harmonic map
		 *	\param epsilon error threshould
		 */	
		void _iterative_map( double threshould = 5e-4 );

	protected:
		/*!	fix the boundary vertices to the unit circle
		 *  using arc length parameter
		 */
		void _set_boundary();
		/*!	The input surface mesh
		 */
		M* m_pMesh;
		/*!	The boundary of m_pMesh
		 */
		typename M::CBoundary m_boundary;
		
	};


/*!	CHarmonicMapper constructor 
*	Count the number of interior vertices, boundary vertices and the edge weight
*
*/
template<typename M>
CDiskHarmonicMapper<M>::CDiskHarmonicMapper( M* pMesh ): m_pMesh( pMesh ), m_boundary( m_pMesh )
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * pV = *viter;
		pV->fixed() = pV->boundary();
	}
	
	//Compute cotangent edge weight
	COperator<M> pC( m_pMesh );
	pC._embedding_2_metric();  //convert embedding to metric
	pC._metric_2_angle();	   //convert metric to angle
	pC._angle_2_Laplace();	   //convert angle to cotangent edge weight

}

//Destructor
/*!
 *	CHarmonicMapper destructor
 */
template<typename M>
CDiskHarmonicMapper<M>::~CDiskHarmonicMapper()
{
}


//Set the boundary vertices to the unit circle
/*!
 *	Fix the boundary using arc length parameter
 */
template<typename M>
void CDiskHarmonicMapper<M>::_set_boundary()
{
	//get the boundary half edge loop
	std::vector<M::CLoop*> & pLs =  m_boundary.loops();
	assert( pLs.size() == 1 );
	M::CLoop * pL = pLs[0];
	std::list<M::CHalfEdge*> & pHs = pL->halfedges();
	
	//compute the total length of the boundary
	double sum = 0;
	for( std::list<M::CHalfEdge*>::iterator hiter = pHs.begin(); hiter != pHs.end(); hiter ++ )
	{
		M::CHalfEdge * ph = *hiter;
		M::CEdge * pE = m_pMesh->halfedgeEdge( ph );
		sum += pE->length();
	}

	//parameterize the boundary using arc length parameter
	double l = 0;
	for( std::list<M::CHalfEdge*>::iterator hiter = pHs.begin(); hiter != pHs.end(); hiter ++ )
	{
		M::CHalfEdge * ph = *hiter;
		M::CEdge * pE = m_pMesh->halfedgeEdge( ph );
		l += pE->length();
		double ang = l/sum * 2.0 * PI;
		M::CVertex * pV = m_pMesh->halfedgeTarget( ph );
		pV->huv()= CPoint2( (cos(ang)+1.0)/2.0, (sin(ang)+1.0)/2.0 );
	}
}

//Compute the harmonic map with the boundary condition, direct method
/*!	Compute harmonic map using direct method
*/
template <typename M>
void CDiskHarmonicMapper<M>::_map()
{
	//fix the boundary
	_set_boundary();

	for( int k = 0; k < 2; k ++ )
	{
		CLaplace<M> L( m_pMesh );

		for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
		{
			M::CVertex * pV = *viter;
			if( !pV->fixed() ) continue;
			pV->u() = pV->huv()[k];
		}
		
		L.solve();

		for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
		{
			M::CVertex * pV = *viter;
			if( pV->fixed() ) continue;
			pV->huv()[k] = pV->u();
		}
	}
}

//Compute the harmonic map with the boundary condition, iterative method
/*!	Iterative method compute harmonic map
*	\param epsilon error threshould
*/
template<typename M>
void CDiskHarmonicMapper<M>::_iterative_map( double epsilon )
{
	//fix the boundary
	_set_boundary();

	//move interior each vertex to its center of neighbors
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * pV = *viter;
		if( pV->fixed() ) continue;
		
		pV->huv() = CPoint2(0,0);
	}
	
	while( true )
	{
		double error = -1e+10;
		//move interior each vertex to its center of neighbors
		for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
		{
			M::CVertex * pV = *viter;
			if( pV->fixed() ) continue;
			
			double  sw = 0;
			CPoint2 suv(0,0);
			for( M::VertexVertexIterator vviter(pV); !vviter.end(); vviter ++ )
			{
				M::CVertex * pW = *vviter;
				M::CEdge   * pE = m_pMesh->vertexEdge( pV, pW );
				double w = pE->weight();
				sw += w;
				suv = suv + pW->huv() * w;
			}
			suv /= sw;
			double verror = (pV->huv()-suv).norm();
			error = (verror > error )?verror:error; 
			pV->huv() = suv;
		}
#ifdef _HARMONIC_MAP_DEBUG_
		printf("Current max error is %f\n", error );
#endif
		if( error < epsilon ) break;
	}	
}




} //namespace Holomorphy
} //namespace MeshLib

#endif

