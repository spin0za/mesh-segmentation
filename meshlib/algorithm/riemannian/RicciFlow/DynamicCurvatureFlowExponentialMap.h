/*! \file DynamicCurvatureFlowExponentialMap.h
 * \brief Dynamic Curvature Flow Exponential Map
 *  \author David Gu
 *  \date   documented on 10/13/2013
 *
 *	Algorithm for Dynamic Curvature Flow Exponential Map
 */

#ifndef _DYNAMIC_CURVATURE_FLOW_EXPONENTIAL_MAP_H_
#define _DYNAMIC_CURVATURE_FLOW_EXPONENTIAL_MAP_H_

#include <map>
#include <vector>


namespace MeshLib
{

namespace RicciFlow
{
	/*! \brief CDynamicCurvatureFlowExponentialMap class
	 *  
	 *  Algorithm for Exponential Map for Dynamic Curvature Flow Embedding result
	 */
	 
 template<typename M>
 class CDynamicCurvatureFlowExponentialMap
  {
  public:
	  /*! \brief CDynamicCurvatureFlowExponentialMap constructor
	   *
	   *  call base class constructor, pCMesh the original mesh, pOMesh the sliced mesh 
	   */
	  CDynamicCurvatureFlowExponentialMap( M * pCMesh, M * pOMesh );

  protected:
	 M * m_pCMesh;
	 M * m_pOMesh;
  };

template<typename M>
CDynamicCurvatureFlowExponentialMap<M>::CDynamicCurvatureFlowExponentialMap( M * pCMesh, M * pOMesh  )
{
	m_pCMesh = pCMesh;
	m_pOMesh = pOMesh;

	for(  M::MeshVertexIterator viter( m_pCMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		pV->valence() = 0;
	}

	for( M::MeshVertexIterator viter( m_pOMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		int id = pV->father();
		M::CVertex * pW = m_pCMesh->idVertex( id );
		pW->valence() ++;
	}

	for( M::MeshVertexIterator viter( m_pOMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		int id = pV->father();
		M::CVertex * pW = m_pCMesh->idVertex( id );
		if( pW->boundary() )
		{
			pV->valence() = pW->valence();
		}
		else
		{
			pV->valence() = 0;
		}
	}


	std::vector<M::CVertex *> corners;
	std::vector<M::CVertex *> anchors;

	M::CBoundary bnd( m_pCMesh );
	int id = 0;
	for( std::vector<M::CLoop* >::iterator liter = bnd.loops().begin(); 
		liter != bnd.loops().end(); liter ++ )
	{
		if( id ++ > 1 ) continue;
		M::CLoop * pL = *liter;

		std::list<M::CHalfEdge*> & pEs = pL->halfedges();
		
		for( std::list<M::CHalfEdge*>::iterator hiter = pEs.begin(); hiter != pEs.end(); hiter ++ )
		{
			M::CHalfEdge * pH = *hiter;
			M::CVertex   * pV = m_pCMesh->halfedgeTarget( pH );
			if( pV->valence() == 2 )
			{
				anchors.push_back( pV );
				break;
			}
		}
	}

	M::CBoundary open_bnd( m_pOMesh );

	for( std::vector<M::CLoop* >::iterator liter = open_bnd.loops().begin();
		liter != open_bnd.loops().end(); liter ++ )
	{
		M::CLoop * open_pL = *liter;

		std::list<M::CHalfEdge*> & open_pEs = open_pL->halfedges();
	
		for( std::list<M::CHalfEdge*>::iterator hiter = open_pEs.begin(); 
			hiter != open_pEs.end(); hiter ++ )
		{
			M::CHalfEdge  * pH = *hiter;
			M::CVertex    * pV = m_pOMesh->halfedgeTarget( pH );
			int id = pV->father();
			M::CVertex    * pW = m_pCMesh->idVertex( id );
			if( pW == anchors[0] )
			{
				corners.push_back( pV );
			}
		}

	}

	for( std::vector<M::CLoop*>::iterator liter = open_bnd.loops().begin();
		liter != open_bnd.loops().end(); liter ++ )
	{
		M::CLoop * open_pL = *liter;

		std::list<M::CHalfEdge*> & open_pEs = open_pL->halfedges();
	
		for( std::list<M::CHalfEdge*>::iterator hiter = open_pEs.begin(); 
			hiter != open_pEs.end(); hiter ++ )
		{
			M::CHalfEdge   * pH = *hiter;
			M::CVertex     * pV = m_pOMesh->halfedgeTarget( pH );
			int id = pV->father();
			M::CVertex     * pW = m_pCMesh->idVertex( id );
			if( pW == anchors[1] )
			{
				corners.push_back( pV );
			}
		}
	}




	 CPoint2 dx = corners[1]->huv() - corners[0]->huv();
	 
	 dx = dx/sqrt( dx[0] * dx[0] + dx[1] * dx[1] );
	 CPoint2 dy( -dx[1], dx[0] );

	 CPoint2 origin = corners[0]->huv();
	
	 for(  M::MeshVertexIterator viter( m_pOMesh ); !viter.end(); viter ++ )
	 {
		 M::CVertex * pV = *viter;
		 CPoint2 uv = pV->huv( );
		 uv = uv - origin;
		 CPoint2 nuv( uv * dx, uv * dy );
		 pV->huv() = nuv;
	 }




	 double width  = (corners[1]->huv() - corners[0]->huv())[0];
	 double height = (corners[2]->huv() - corners[0]->huv())[1];
	
	 width  = fabs( width  );
	 height = fabs( height );


	 for(  M::MeshVertexIterator viter( m_pOMesh ); !viter.end(); viter ++ )
	 {
		 M::CVertex * pV = *viter;
		 CPoint2 p = pV->huv();

		 p = p/width * 2 *PI;
			
		 double rad = exp( height )/exp( p[1] );
		 double ang = p[0];

		 CPoint2 uv( rad * cos(ang), rad * sin(ang ));
		 pV->huv() = uv;
	}

	for(  M::MeshVertexIterator viter( m_pOMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		int id = pV->father();
		M::CVertex * pW = m_pCMesh->idVertex( id );
		pW->huv() = pV->huv();
	}

	if( mag2( anchors[0]->huv( ) ) > mag2( anchors[1]->huv()) )
	{
			return;
	}

	for(  M::MeshVertexIterator viter( m_pOMesh ); !viter.end(); viter ++ )
	 {
		 M::CVertex * pV = *viter;
		 CPoint2 p = pV->huv();

		 double r = mag( p );			
		 CPoint2 uv( p[0]/(r*r), p[1]/(r*r) );
		 pV->huv() = uv;
	}


	for(  M::MeshVertexIterator viter( m_pOMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		int id = pV->father();
		M::CVertex * pW = m_pCMesh->idVertex( id );
		pW->huv() = pV->huv();
	}

	//write_uv( mesh );
	//mesh.write_m( argv[3] );

	//return 0;
};


} //namespace RicciFlow
} //namespace MeshLib

#endif  