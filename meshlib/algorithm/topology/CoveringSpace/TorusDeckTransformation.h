/*!
*      \file TorusDeckTransformation.h
*      \brief Compute Deck Transformation Group
*	   \author David Gu
*      \date Documented 03/31/2011
*
*/

#ifndef  _TORUS_DECK_TRANSFORMATION_H_
#define  _TORUS_DECK_TRANSFORMATION_H_

#include <map>
#include <queue>
#include <vector>
#include <algorithm>

#include "Riemannian/RicciFlow/RicciFlowMesh.h"
#include "DeckTransformation.h"

namespace MeshLib
{

namespace Topology
{

/*! Compuare two Ricci Vertices according to their father field */
template<typename V>
class CompareRicciFlowVertex 
{
public:
    bool operator()(V * pV1, V * pV2)
    {
       return (pV1->father() > pV2->father() );
    }
};

/*!
*	\brief CTorusDeckTransformation
*	
*	Deck Transformation Group for Torus
* 
*/
template<typename M>
class CTorusDeckTransformation: public CDeckTransformation<M>
{

public:
	/*!
	 *	CTorusDeckTransformation constructor
	 *  \param pMesh    closed mesh
	 *  \param pDomain  foundamental domain
	 */
	CTorusDeckTransformation( M * pMesh, M * pDomain ): CDeckTransformation( pMesh, pDomain) {};
	/*!
	 *	CTorusDeckTransformation destructor
	 */
	~CTorusDeckTransformation() {};

	void _calculate_generators();

	std::vector<CPoint2> & generators() { return m_translations; };

protected:

	bool _locate_corners( std::vector<typename M::CVertex*> & corners );
	
	std::vector<CPoint2> m_translations;

};

template<typename M>
inline bool CTorusDeckTransformation<M>::_locate_corners( std::vector<typename M::CVertex*> & corners )
{

	std::priority_queue<M::CVertex*, std::vector<M::CVertex*>, CompareRicciFlowVertex<M::CVertex> > vqueue;

	for( unsigned int i = 0; i < m_boundary.loops().size(); i ++ )
	{
		M::CLoop * pL = m_boundary.loops()[i];

		for( std::list<M::CHalfEdge*>::iterator hiter = pL->halfedges().begin(); hiter != pL->halfedges().end(); hiter ++ )
		{
			M::CHalfEdge * pH = *hiter;
			M::CVertex * pV = m_pDomain->halfedgeTarget( pH );
			vqueue.push( pV );
		}
	
		std::queue<M::CVertex*> corner_queue;

		for( int j = 0; j < 3 && !vqueue.empty(); j ++ )
		{
			M::CVertex *pV = vqueue.top();
			vqueue.pop();
			corner_queue.push( pV );
		}
		
		if( corner_queue.size() < 3 ) continue; //the boundary only has 3 vertices

		while( !vqueue.empty() )
		{
			M::CVertex * pV = vqueue.top();
			vqueue.pop();
			corner_queue.push( pV );
						
			bool beCorners = true;
			
			for( int j = 0; j < 4; j ++ )
			{
				M::CVertex * pV = corner_queue.front();
				corner_queue.pop();
				corner_queue.push( pV );
				if( pV->father() != corner_queue.front()->father() ) beCorners = false;
			}

			if( beCorners )
			{
				printf("Locate corners\n");
				while( !corner_queue.empty() )
				{
					M::CVertex * pV = corner_queue.front();
					corner_queue.pop();
					corners.push_back( pV );
					printf("Vertex %d ->Father is %d\n", pV->id(), pV->father());
				}

				return true;
			}

			corner_queue.pop();
		}
	}

	return false;
};

template<typename M>
inline void CTorusDeckTransformation<M>::_calculate_generators()
{
	std::vector<M::CVertex*> corners;

	if( !_locate_corners( corners ) )	
	{
		std::cout <<"Error: Could not find all the four corners." << std::endl;
		return;
	}

	std::vector<M::CVertex*> vsequence;

	for( unsigned int i = 0; i < m_boundary.loops().size(); i ++ )
	{
		M::CLoop * pL = m_boundary.loops()[i];
		
		bool interior_boundary = true;
		for( std::list<M::CHalfEdge*>::iterator hiter = pL->halfedges().begin(); hiter != pL->halfedges().end(); hiter ++ )
		{
			M::CHalfEdge * pH = *hiter;
			M::CVertex * pV = m_pDomain->halfedgeTarget( pH );
			if( pV == corners[0] ) 
			{
				interior_boundary = false;
				break;
			}
		}
		if( interior_boundary ) continue;

		//this is the boundary loop
		for( std::list<M::CHalfEdge*>::iterator hiter = pL->halfedges().begin(); hiter != pL->halfedges().end(); hiter ++ )
		{
			M::CHalfEdge * pH = *hiter;
			M::CVertex * pV = m_pDomain->halfedgeTarget( pH );
			vsequence.push_back( pV );
		}
	}

	while( true )
	{
		
		M::CVertex * pV = vsequence.front();
		std::vector<M::CVertex*>::iterator  viter;
		viter = std::find ( corners.begin(), corners.end(), pV );	

		if( viter != corners.end() ) break;

		vsequence.erase( vsequence.begin() );
		vsequence.push_back( pV );
	}

	std::vector<M::CVertex*> sequence_corners;

	for( unsigned int i = 0; i < vsequence.size(); i ++ )
	{
		M::CVertex * pV = vsequence[i];
		std::vector<M::CVertex*>::iterator  viter;
		viter = std::find ( corners.begin(), corners.end(), pV );
		if( viter == corners.end() ) continue;
		sequence_corners.push_back( pV );		
	}

	std::cout << "Found "<< corners.size() << " Sequencial Corners" << std::endl;

	for( unsigned int i = 0; i < 4; i ++ )
	{
		CPoint2 uv = sequence_corners[i]->huv();
		std::cout << "Father is "<< corners.size() << " Sequencial Corners" << std::endl;
	}

	CPoint2  Ta = sequence_corners[1]->huv() - sequence_corners[0]->huv();
	CPoint2  Tb = sequence_corners[2]->huv() - sequence_corners[1]->huv();
	
	m_translations.push_back( Ta );
	m_translations.push_back( Tb );

	for( unsigned int i = 0; i < m_translations.size(); i ++ )
	{
		CPoint2 T = m_translations[i];
		std::cout <<  T[0] << " " << T[1] << std::endl;
	}

};

} //namespace Topology
} //namespace MeshLib

#endif  _TORUS_DECK_TRANSFORMATION_H_