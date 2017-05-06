/*! \file CanonicalBasis.h
*   \brief Algorithm for computing the canonical fundamental group generators
*	\author David Gu
*   \date Documented on 04/07/2013
*/

#ifndef _CANONICAL_BASIS_H_
#define _CANONICAL_BASIS_H_

#include <queue>
#include <deque>
#include <set>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream
#include <fstream>
#include "Mesh/boundary.h"
#include "CanonicalBasisMesh.h"

namespace MeshLib
{
namespace Topology
{

/*! Compuare two Dijkstra Vertices according to their father traits */

class CompareEdgeAncestor 
{
public:
    bool operator()(CCanonicalBasisEdge * pE1, CCanonicalBasisEdge * pE2)
    {
		CCBMesh * pM = (CCBMesh * ) CCanonicalBasisEdge::m_pMesh;

		CCanonicalBasisVertex * pV0 = pM->edgeVertex1( pE1 );
		CCanonicalBasisVertex * pV1 = pM->edgeVertex2( pE1 );

		int id0 = (pV0->ancestor() < pV1->ancestor())?pV0->ancestor() : pV1->ancestor();
		int id1 = (pV0->ancestor() > pV1->ancestor())?pV0->ancestor() : pV1->ancestor();

		CCanonicalBasisVertex * pW0 = pM->edgeVertex1( pE2 );
		CCanonicalBasisVertex * pW1 = pM->edgeVertex2( pE2 );

		int iD0 = (pW0->ancestor() < pW1->ancestor())?pW0->ancestor() : pW1->ancestor();
		int iD1 = (pW0->ancestor() > pW1->ancestor())?pW0->ancestor() : pW1->ancestor();

		if( id0 > iD0 ) return true;
		if( id0 < iD0 ) return false;

		return (id1 > iD1);
	}
};

	
class CCanonicalBasisSegment
{
public:
	CCanonicalBasisSegment( CCBMesh * pMesh )
	{
		m_pMesh = pMesh;
	};
	std::vector<CHalfEdge*> & hes()
	{
		return m_hes;
	}
	CCanonicalBasisSegment* & dual() 
	{
		return m_dual;
	}
	bool angle_ascending();

	CCanonicalBasisVertex * start()
	{
		if( m_hes.empty() ) return NULL;

		CHalfEdge * pH = m_hes[0];
		CCanonicalBasisVertex * pV = m_pMesh->halfedgeSource( pH );
		return pV;
	};

	//output to a txt file
	void write( const char * name );
protected:
	CCBMesh * m_pMesh;
	std::vector<CHalfEdge*> m_hes;
	CCanonicalBasisSegment *m_dual;
};

inline bool CCanonicalBasisSegment::angle_ascending()
{
	CCanonicalBasisVertex * pB = m_pMesh->halfedgeSource( m_hes[0] );
	std::vector<double> areas;

	for( size_t i = 0; i < m_hes.size(); i ++ )
	{
		CHalfEdge * pH = m_hes[i];
		CCanonicalBasisVertex * pV = m_pMesh->halfedgeSource( pH );
		CCanonicalBasisVertex * pW = m_pMesh->halfedgeTarget( pH );
		double area = (pV->point()^pW->point())* CPoint(0,0,1);
		areas.push_back( area );
	}
	
	for( size_t i = 0; i < areas.size(); i ++ )
	{
		if( areas[i] < -1e-8 ) return false;
	}

	return true;
};

inline void CCanonicalBasisSegment::write( const char * name )
{
	std::ofstream fp;
	fp.open(name);

	for( size_t i = 0; i < m_hes.size(); i ++ )
	{
		CHalfEdge * pH = m_hes[i];
		CCanonicalBasisVertex * pV = m_pMesh->halfedgeSource( pH );
		CCanonicalBasisVertex * pW = m_pMesh->halfedgeTarget( pH );
		fp << pV->ancestor() << " " << pW->ancestor() << std::endl;
	}
	fp.close();
}

class CHandle
{
public:
	CHandle( CCBMesh * pMesh ):m_pMesh(pMesh){};
	~CHandle()
	{ 
		for( size_t i = 0; i < m_segments.size(); i ++ )
		{
			delete m_segments[i];
			m_segments.clear();
		}
	};
	std::set<CCanonicalBasisVertex*> & anchors() { return m_anchors; };
	std::vector<CCanonicalBasisSegment*> & segments() { return m_segments; };

	void sort();

	CCanonicalBasisVertex * start()
	{
		if( m_segments.empty() ) return NULL;
		return m_segments[0]->start();
	}
protected:
	CCBMesh * m_pMesh;
	std::vector<CCanonicalBasisSegment*> m_segments;
	std::set<CCanonicalBasisVertex*> m_anchors;
};

//sort the segments in this handle
//in the first segment, the argument is monotonously ascending
//there are 4 segments, one ascending; one descending; other two turning at the corners

inline void CHandle::sort()
{
	std::deque<CCanonicalBasisSegment*> segs;
	for( size_t i = 0; i < m_segments.size(); i ++ )
	{
		segs.push_back( m_segments[i] );
	}
	//find the first segment, 
	while( true )
	{
		CCanonicalBasisSegment * pS = segs.front();
		if( pS->angle_ascending() ) break;
		segs.pop_front();
		segs.push_back( pS );
	}

	m_segments.clear();

	while( !segs.empty() )
	{
		m_segments.push_back( segs.front() );
		segs.pop_front();
	}
};

/*! Compuare two Handle by the argument of their starting vertex */

class CompareHandle 
{
public:
    bool operator()(CHandle * pH1, CHandle * pH2)
    {

		CCanonicalBasisVertex * pV1 = pH1->start();
		CCanonicalBasisVertex * pV2 = pH2->start();

		double ang1 = atan2( pV1->point()[1], pV1->point()[0] );
		double ang2 = atan2( pV2->point()[1], pV2->point()[0] );

		return (ang1 > ang2);
	}
};



/*! \brief CCanonicalBasis class
*
*	Compute the canonical fundamental group basis
*
*/ 
class CCanonicalBasis
{
public:
	/*! CDijkstra constructor
	* \param pMesh input closed mesh 
	* \param pDomain open mesh, pMesh sliced along a loop
	*/
	CCanonicalBasis( CCBMesh * pDomain ):m_boundary( pDomain )
	{
		m_pMesh   = pDomain;
		CCanonicalBasisEdge::m_pMesh = (void*) m_pMesh;

	};

	/*! CDijkstra destructor */
	~CCanonicalBasis()
	{
		for( size_t i = 0; i < m_handles.size(); i ++ )
		{
			delete m_handles[i];
		}
		m_handles.clear();
	};

	/*! Shortest Basis */
	void _shortest_basis( CDijkstraVertex * root );

	/*! 3. \beta_k's are a_k's, then compute the conjugate b_k */
	int shortest_path( const char * loop_file, const char * father_loop_file );

	void _canonical_homotopy_basis();
private:

	/*! Closed mesh */
	CCBMesh * m_pMesh;
	/*! boundary of the open mesh */
	CCBMesh::CBoundary m_boundary;

protected:
	void _construct_handle( CCBMesh::CLoop * pL );
	std::vector<CHandle*> m_handles;
};


inline void CCanonicalBasis::_construct_handle( CCBMesh::CLoop * pL )
{
	CHandle * pHandle = new CHandle( m_pMesh );
	assert( pHandle );

	m_handles.push_back( pHandle );

	std::priority_queue<CCanonicalBasisVertex*, std::vector<CCanonicalBasisVertex*>, MeshLib::Topology::CompareVertexAncestor>  vqueue;

	//find the four anchors

	for( std::list<CHalfEdge*>::iterator hiter = pL->halfedges().begin(); hiter != pL->halfedges().end(); hiter ++ )
	{
		CHalfEdge * pH = *hiter;
		CCanonicalBasisVertex * pT = m_pMesh->halfedgeTarget(pH);
		vqueue.push( pT );
	}

	while( !vqueue.empty() )
	{
		CCanonicalBasisVertex * pV = vqueue.top();
		vqueue.pop();
		CCanonicalBasisVertex * pW = vqueue.top();
		vqueue.pop();

		if( !vqueue.empty() && pW->ancestor() == vqueue.top()->ancestor() )
		{
			std::cout << "Anchor" << std::endl;
			
			pHandle->anchors().insert( pV );
			pHandle->anchors().insert( pW );

			for( int i = 0; i < 2; i ++ )
			{
				pV = vqueue.top();
				vqueue.pop();
				pHandle->anchors().insert( pV );
			}
		}
	}

	std::cout << "There are " << pHandle->anchors().size() << " anchors" << std::endl;

	std::deque<CHalfEdge*> dq;

	for( std::list<CHalfEdge*>::iterator hiter = pL->halfedges().begin(); hiter != pL->halfedges().end(); hiter ++ )
	{
		CHalfEdge * pH = *hiter;
		dq.push_back( pH );
	}

	//find the first half-edge
	while( true )
	{
		CHalfEdge * pH = dq.front();
		CCanonicalBasisVertex * pV = m_pMesh->halfedgeSource( pH );
		if( pHandle->anchors().find( pV ) != pHandle->anchors().end() ) break;
		dq.pop_front();
		dq.push_back( pH );			
	}


	//split the boundary loop to 4 segments
	while( !dq.empty() )
	{
		CCanonicalBasisSegment * pS = new CCanonicalBasisSegment( m_pMesh );
		assert( pS != NULL );
		pHandle->segments().push_back( pS );

		//find the tail
		while( true && !dq.empty() )
		{
			CHalfEdge * pH = dq.front();
			dq.pop_front();

			pS->hes().push_back( pH );	
			
			CCanonicalBasisVertex * pW = m_pMesh->halfedgeTarget( pH );
			if( pHandle->anchors().find( pW ) != pHandle->anchors().end() ) break;
		}
	}

	pHandle->segments()[0]->dual() = pHandle->segments()[2];
	pHandle->segments()[2]->dual() = pHandle->segments()[0];
	pHandle->segments()[1]->dual() = pHandle->segments()[3];
	pHandle->segments()[3]->dual() = pHandle->segments()[1];

	pHandle->sort();
};

inline void CCanonicalBasis::_canonical_homotopy_basis()
{
	std::vector<CCanonicalBasisSegment*> segments;

	std::priority_queue<CCanonicalBasisVertex*, std::vector<CCanonicalBasisVertex*>, MeshLib::Topology::CompareVertexAncestor>  vqueue;

	for( unsigned int i = 0; i < m_boundary.loops().size(); i ++ )
	{
		CCBMesh::CLoop * pL = m_boundary.loops()[i];
		_construct_handle( pL );
	}

	//sort the handles, using slit map coordinates
	std::priority_queue<CHandle*, std::vector<CHandle*>, MeshLib::Topology::CompareHandle>  hqueue;
	
	for( size_t i = 0; i < m_handles.size(); i ++)
	{
		CHandle * pH = m_handles[i];
		hqueue.push( pH );
	}

	m_handles.clear();

	while( ! hqueue.empty() )
	{
		CHandle * pH = hqueue.top();
		hqueue.pop();
		m_handles.push_back( pH );
	}


/*
	std::priority_queue<CCanonicalBasisEdge*, std::vector<CCanonicalBasisEdge*>, MeshLib::CompareEdgeAncestor>  equeue;
	
	for( unsigned int i = 0; i < m_boundary.loops().size(); i ++ )
	{
		CCBMesh::CLoop * pL = m_boundary.loops()[i];
		
		for( std::list<CHalfEdge*>::iterator hiter = pL->halfedges().begin(); hiter != pL->halfedges().end(); hiter ++ )
		{
			CHalfEdge * pH = *hiter;
			CCanonicalBasisEdge * pE = m_pMesh->halfedgeEdge( pH );
			equeue.push( pE );
		}

		while( !equeue.empty() )
		{
			CCanonicalBasisEdge * pE0 = equeue.top();
			equeue.pop();
			CCanonicalBasisEdge * pE1 = equeue.top();
			equeue.pop();
			
			std::cout << "(" << m_pMesh->edgeVertex1( pE0 )->ancestor() <<"," << m_pMesh->edgeVertex2( pE0 )->ancestor() << ")" << std::endl;
			std::cout << "(" << m_pMesh->edgeVertex1( pE1 )->ancestor() <<"," << m_pMesh->edgeVertex2( pE1 )->ancestor() << ")" << std::endl;
		}
	}
*/

	//output

	for( size_t  i = 0; i < m_handles.size(); i ++ )
	{
		std::stringstream ss;
		ss << "a_" << i << ".txt";
		m_handles[i]->segments()[0]->write( ss.str().c_str() ); 

		std::stringstream ns;
		ns << "b_" << i << ".txt";
		m_handles[i]->segments()[1]->write( ns.str().c_str() );
	}
};


} //namespace Topology
} //namespace MeshLib
#endif