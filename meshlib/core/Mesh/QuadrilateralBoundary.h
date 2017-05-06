/*!
*      \file QuadrilateralBoundary.h
*      \brief Trace boundary loops
*	   \author David Gu
*      \date 11/23/2012
*
*/
//the four corners of the boundary should be marked as ``marker''

#ifndef _QUADRILATERAL_BOUNDARY_H_
#define _QUADRILATERAL_BOUNDARY_H_

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <list>
#include <set>

#include "../Mesh/BaseMesh.h"
#include "../Mesh/iterators.h"
#include "boundary.h"

namespace MeshLib
{
/*!
	\brief CLoop Boundary loop  class.	
	\tparam CVertex Vertex type
	\tparam CEdge   Edge   type
	\tparam CFace   Face   type
	\tparam CHalfEdge HalfEdge type
*/

template<typename M, typename V, typename E, typename F, typename H>
class CQuadrilateralBoundary
{
public:
	/*!
		Constructor of the CQuadrilateralBoundary
		\param pMesh  pointer to the current mesh
	*/
	 CQuadrilateralBoundary( M * pMesh );
	 /*!
		Destructor of CQuadrilateralBoundary.
	 */
	~CQuadrilateralBoundary();
	/*!	
	 * get array of segments	
	 */
	std::vector<std::list<H*>*> & segments() { return m_segments; };

protected:
	/*!
	 *	Pointer to the current mesh.
	 */
	M		* m_pMesh;
	/*!
	 *	Segments
	 */
	std::vector<std::list<H*>*> m_segments;
};

/*!
	CBoundary constructor
	\param pMesh the current mesh
*/
template<typename M, typename V, typename E, typename F, typename H>
CQuadrilateralBoundary<M,V,E,F,H>::CQuadrilateralBoundary( M * pMesh )
{
	m_pMesh = pMesh;
	CBoundary<V,E,F,H> boundary( pMesh );

	//get the boundary half edge loop
	std::vector<M::CLoop*> & pLs =  boundary.loops();
	assert( pLs.size() == 1 );
	M::CLoop * pL = pLs[0];
	std::list<H*>   pHs; 
	for( std::list<H*>::iterator hiter = pL->halfedges().begin(); hiter != pL->halfedges().end(); hiter ++ )
	{
		H * ph = *hiter;
		pHs.push_back( ph );
	}
	
	std::list<V*> corners;
	for( std::list<H*>::iterator hiter = pHs.begin(); hiter != pHs.end(); hiter ++ )
	{
		H * ph = *hiter;
		V * pv = m_pMesh->halfedgeSource( ph );
		//int valence = 0;
		//for( M::VertexEdgeIterator veiter( pv ); !veiter.end(); veiter ++ )
		//{
		//	E * pe = *veiter;
		//	valence += ( pe->sharp() )?1:0;
		//}
		//if( valence > 2 ) corners.push_back( pv );
		if( pv->isMarker() ) corners.push_back( pv );
	}
	
	std::cout <<"CQuadrilateralBoundary::CQuadrilateralBoundary find "<< corners.size()<< std::endl;
	
	if( corners.size() != 4 ) 
	{
		std::cerr << "There should be 4 corners" << std::endl;
		return;
	}

	//find the lower left corner
	V * pLL = corners.front();

	for( std::list<V*>::iterator viter = corners.begin(); viter != corners.end(); viter ++ )
	{
		V * pv = *viter;
		if( pv->z().real() < pLL->z().real() )
		{
			pLL = pv;
			continue;
		}
		if( pv->z().real() > pLL->z().real() ) continue;

		if( pv->z().imag() < pLL->z().imag() )
			pLL = pv;
	}

	while( true )
	{
		H * ph = pHs.front();
		V * pv = m_pMesh->halfedgeSource( ph );
		if( pv == pLL ) break;
		pHs.pop_front();
		pHs.push_back( ph );
	}

	std::list<H*> * pSeg = NULL;
	while( !pHs.empty() )
	{
		H * ph = pHs.front();
		V * pv = m_pMesh->halfedgeSource( ph );
		pHs.pop_front();

		if( std::find( corners.begin(), corners.end(), pv ) != corners.end() ) 
		{
			pSeg = new std::list<H*>;
			assert( pSeg != NULL );
			m_segments.push_back( pSeg );
		}
		pSeg->push_back( ph );
	}
};

/*!
	CQuadrilateralBoundary destructor
*/
template<typename M, typename V, typename E, typename F, typename H>
CQuadrilateralBoundary<M,V,E,F,H>::~CQuadrilateralBoundary()
{
	for( size_t i = 0; i < m_segments.size(); i ++ )
	{
		std::list<H*> * pSeg = m_segments[i];
		delete pSeg;
	}
	m_segments.clear();
};

}
#endif

