/*!
*      \file Cohomology.h
*      \brief Algorithm for base 1-forms for cohomology group (for closed surfaces)
*	   \author David Gu
*      \date Document 12/26/2010
*
*		Computing base 1-forms for cohomology group(for closed surfaces)
*
*/

/********************************************************************************************************************************
*
*      Cohomology Form Class
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       Cohomology Group base 1-Form
* 
*       David Gu December 26, 2010,  gu@cs.stonybrook.edu
*
*
*      Input:
*         
*           Original mesh, the mesh cut through a homology basis
*
*      Output:
*
*           Closed non-exact closed 1-form. and the mesh with UV coordinates.
*
*********************************************************************************************************************************/


#ifndef _COHOMOLOGY_H_
#define _COHOMOLOGY_H_

#include "CohomologyMesh.h"

namespace MeshLib
{

namespace Topology
{
	/*! \brief CCohomology class
	*  
	*  Algorithm for computing closed 1-form basis for the 1st cohomology group
	*/
  template<typename M>
  class CCohomology
  {
  public:
	  /*! CCohomology constructor
	  *
	  *	\param pMesh input closed mesh
	  * \param pWMesh sliced open mesh
	  */
	
    CCohomology( M * pMesh, M * pWMesh );
	/*! CCohomology destructor
	*/

    ~CCohomology();

	/*!	Compute harmonic closed forms on the input mesh
	*/
    void calculate_closed_form();

  protected:
	/*! Input closed mesh */	
    M * m_pMesh;
	/*! inpute mesh, sliced open along the shortest paths connecting two boundary loops */    
	M * m_pWMesh;

	typename M::CBoundary m_boundary;

	/*! Compute the closed 1-form*/
    void _closed_1_form();
	/*! Integrate the closed form on the mesh */
    void _integrate();

  };

template<typename M>
CCohomology<M>::CCohomology( M * pMesh, M * pWMesh ):m_boundary( pWMesh )
{
  m_pMesh  = pMesh;
  m_pWMesh = pWMesh;
};

template<typename M>
CCohomology<M>::~CCohomology()
{
};

template<typename M>
void CCohomology<M>::_integrate()
{
	for( M::MeshEdgeIterator eiter( m_pWMesh ); !eiter.end(); eiter ++ )
  {
	  M::CEdge * we = *eiter;

	  M::CVertex * w1 = m_pWMesh->edgeVertex1( we );
	  M::CVertex * w2 = m_pWMesh->edgeVertex2( we );

      int id1 = w1->father();
      int id2 = w2->father();

	  M::CVertex * v1 = m_pMesh->idVertex( id1 );
	  M::CVertex * v2 = m_pMesh->idVertex( id2 );

	  M::CEdge * e = m_pMesh->vertexEdge( v1, v2 );

      if( m_pMesh->edgeVertex1( e ) == v1 )
      {
          we->du() = e->du();
      }
      else
      {  
        we->du() = - e->du();
      }
  }

	M::CVertex * head = NULL;
  for( M::MeshVertexIterator viter( m_pWMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
      v->touched() = false;
      head = v;
  }

  std::queue<M::CVertex*> vqueue;
  head->touched() = true;
  vqueue.push( head );
  head->uv() = CPoint2(0,0.27);

  while( !vqueue.empty() )
  {
	  M::CVertex * head = vqueue.front();
    vqueue.pop();

	for( M::VertexEdgeIterator veiter(  head ); !veiter.end(); veiter ++ )
    {
		M::CEdge * e = *veiter;
      double du = e->du();

	  M::CVertex * v1 = m_pWMesh->edgeVertex1( e );
	  M::CVertex * v2 = m_pWMesh->edgeVertex2( e );

      if( v1 == head )
      {
        if( v2->touched() ) continue;
        v2->touched() = true;
        vqueue.push( v2 );
        v2->uv() = v1->uv() + CPoint2( du,0 );
      }
      else
      {
        assert( v2 == head );
        if( v1->touched() ) continue;
        v1->touched() = true;
        vqueue.push( v1 );
        v1->uv() = v2->uv() - CPoint2( du,0 );
      }
    }
  }
  
  double umin = 1e+10;
  double umax = -1e+10;

  for( M::MeshVertexIterator viter( m_pWMesh );  !viter.end(); viter ++ )
  {
	  M::CVertex* v = *viter;
  	  CPoint2 p = v->uv();
      umax = (umax > p[0] )? umax: p[0];
      umin = (umin < p[0] )? umin: p[0];
  }

  for( M::MeshVertexIterator viter( m_pWMesh );  !viter.end(); viter ++ )
  {
	  M::CVertex* v = *viter;
  	  CPoint2 p = v->uv();
      p[0] = ( p[0] - umin )/(umax-umin);
      v->uv() = p;
  }
}

template<typename M>
void CCohomology<M>::calculate_closed_form()
{
    _closed_1_form();
    _integrate();
}

template<typename M>
void CCohomology<M>::_closed_1_form()
{
  for( M::MeshVertexIterator viter( m_pWMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * pV = *viter;
      pV->u() = 0;
  }
  
  std::vector<M::CLoop*> & loops = m_boundary.loops();
  assert( loops.size() == 2 );
  
	M::CLoop * pL = loops[0];
	for( std::list<M::CHalfEdge*>::iterator hiter = pL->halfedges().begin(); hiter != pL->halfedges().end(); hiter ++ )
	{
		M::CHalfEdge	* pH = *hiter;
		M::CVertex      * pV = m_pWMesh->halfedgeTarget( pH );
		pV->u() = 1.0;
	}

  for( M::MeshEdgeIterator eiter( m_pWMesh ); !eiter.end(); eiter ++ )
  {
	  M::CEdge * we = *eiter;
	  M::CVertex * w1 = m_pWMesh->edgeVertex1( we );
	  M::CVertex * w2 = m_pWMesh->edgeVertex2( we );
	  we->du() = w2->u()-w1->u();

    int id1 =  w1->father();
    int id2 =  w2->father();

	M::CVertex * v1 = m_pMesh->idVertex( id1 );
	M::CVertex * v2 = m_pMesh->idVertex( id2 );

	M::CEdge * e = m_pMesh->vertexEdge( v1, v2 );

    if( m_pMesh->edgeVertex1( we ) == v1 )
    {
        e->du() = we->du();
    }
    else
    {
        e->du() = -we->du();
    }
  }
};

} //namespace Topology
} //namespace MeshLib
#endif