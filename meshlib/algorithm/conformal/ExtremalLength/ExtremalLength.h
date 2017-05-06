/*!
*      \file ExtremalLength.h
*      \brief Algorithm for extremal length
*	   \author David Gu
*      \date Document 12/29/2010
*
*		Computing extremal length
*
*/

#ifndef _EXTREMAL_LENGTH_H_
#define _EXTREMAL_LENGTH_H_

#include  <math.h>
#include <queue>
#include "ExtremalLengthMesh.h"
#include "Operator/Operator.h"
#include "Laplace/Laplace.h"

namespace MeshLib
{
namespace Holomorphy
{
	/*! \brief CExtremalLength class
	*  
	*  Algorithm for computing extremal length
	*/
  template<typename M>
  class CExtremalLength
  {
  public:
	  /*! CExtremalLength constructor
	  *
	  *	\param pMesh input mesh, topological disk
	  */
     CExtremalLength( M * pMesh );
	  /*! CExtremalLength destructor
	  */
    ~CExtremalLength();

	/*!	Compute harmonic forms on the input mesh
	* \param prefix for output mesh file	
	*/
    void calculate_harmonic_exact_form( const char * prefix);

  protected:
    /*!	input mesh
	 */
	M * m_pMesh;
	/*! the boundary of the input mesh
	 */
	typename M::CBoundary m_boundary;
	/*!	Segment the boundary to four sides
	 */
	void _segment_boundary();
	/*! Compute the exact harmonic form, the value on top is one, the value on bottom is zero
	 * \param top top side of the boundary sgements
	 * \param bottom bottom side of the boundary segments
	 */
	void _harmonic_exact_form( std::vector<typename M::CVertex*> & top, std::vector<typename M::CVertex*> & bottom);
	/*! number of interior vertices */
	int  m_interior_vertices;
	/*! number of boundary vertices */
	int  m_boundary_vertices;

	/*! four boundary segments of the topological quadralaterial
	*/
	std::vector<std::vector<typename M::CVertex*>*> m_segments;

  };


template<typename M>
CExtremalLength<M>::CExtremalLength( M * pMesh ):m_pMesh( pMesh ), m_boundary( m_pMesh )
{

	_segment_boundary();

	MeshLib::COperator<M> pC( m_pMesh );
	pC._embedding_2_metric();
	pC._metric_2_angle();
	pC._angle_2_Laplace();	//compute edge weight

}

template<typename M>
CExtremalLength<M>::~CExtremalLength()
{
	for( size_t i = 0 ; i < m_segments.size(); i ++ )
	{
		std::vector<M::CVertex*> * pSeg = m_segments[i];
		delete pSeg;
	}
}


template<typename M>
void CExtremalLength<M>::calculate_harmonic_exact_form(const char * prefix )
{

  for( int k = 0; k < 2; k ++ )
  {
		
	  _harmonic_exact_form( *m_segments[k], *m_segments[k+2] );

	 std::stringstream iss;
	 iss << prefix << "_" << k << ".du.m" ;
	 m_pMesh->write_m( iss.str().c_str() );
  }
}

template<typename M>
void CExtremalLength<M>::_harmonic_exact_form( std::vector<typename M::CVertex*>& top_segment,std::vector<typename M::CVertex*>& bottom_segment  )
{

	// initialize vertex index and u value

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * pV = *viter;
		pV->fixed() = false;
		pV->u() = 0.0;
	}


	//set boundary condition

	for( size_t i = 0; i < top_segment.size(); i ++ )
	{
		M::CVertex   * pV = top_segment[i];
		pV->fixed() = true;
		pV->u() = 1.0;
	}

	for( size_t i = 0; i < bottom_segment.size(); i ++ )
	{
		M::CVertex   * pV = bottom_segment[i];
		pV->fixed() = true;
		pV->u() = 0.0;
	}

	CLaplace<M> L( m_pMesh );
	L.solve();

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * pV = *viter;
		if( pV->fixed() )
		{
			pV->uv() = CPoint2( pV->u(), 0.53 );
			continue;
		}
		pV->uv() = CPoint2( pV->u(), 0.53 );
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );
		e->du() = v2->u() - v1->u();
	}
}


template<typename M>
void CExtremalLength<M>::_segment_boundary()
{
	//compute the valence of vertices

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
    {
		M::CVertex * v = *viter;
		v->valence() = 0;
		for( M::VertexEdgeIterator veiter( v ); !veiter.end(); ++ veiter )
		{
			M::CEdge * e = *veiter;
			if( e->sharp() ) v->valence() ++;
		}
		if( v->valence() > 1 && v->boundary() ) 
			  printf("Corner %d\n", v->id() );
	}


	M::CLoop * pL = *( m_boundary.loops().begin() );	
	std::vector<M::CVertex*> vertices;

	for( std::list<M::CHalfEdge*>::iterator hiter = pL->halfedges().begin(); hiter != pL->halfedges().end(); hiter ++ )
	{
		M::CHalfEdge * he = *hiter;
		M::CVertex* pV = m_pMesh->halfedgeTarget( he );
		vertices.push_back( pV );
	}
	
	size_t idx = 0;

	for( idx = 0; idx < vertices.size(); idx ++ )
	{
		M::CVertex * pV = vertices[idx];
		if( pV->valence() > 1 ) break;
	}

	int s = 0;
	for( int i = 0; i < 4; i ++ )
	{
		std::vector<M::CVertex*> * pSeg = new std::vector<M::CVertex*>;
		m_segments.push_back( pSeg );
		assert( pSeg != NULL );

		//starting corner

		M::CVertex * pV = vertices[(idx+s)%vertices.size()];
		pSeg->push_back( pV );

		while( true )
		{
			s++;
			int ind = (idx+s)%vertices.size();
			M::CVertex * pV = vertices[ind];
			pSeg->push_back( pV );
			if( pV->valence() > 1 ) break;
		}

		
	}

	assert( s == vertices.size() );
}


} //namespace Holomorphy
} //namespace MeshLib
#endif _EXTREMAL_LENGTH_H_