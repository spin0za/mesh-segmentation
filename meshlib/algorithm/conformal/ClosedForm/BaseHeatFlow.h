/*!
*      \file BaseHeatFlow.h
*      \brief Algorithm for base harmonic forms (for closed surfaces)
*	   \author David Gu
*      \date Document 12/26/2010
*
*		Computing harmonic forms (for closed surfaces)
*
*/
#ifndef _BASE_HEAT_FLOW_H_
#define _BASE_HEAT_FLOW_H_

#include "HarmonicClosedFormMesh.h"
#include "Operator/Operator.h"
#include "Laplace/Poisson.h"

namespace MeshLib
{

namespace Holomorphy
{

	/*! \brief CBaseHeatFlow class
	*  
	*  Algorithm for computing harmonic forms
	*/
  template<typename M>
  class CBaseHeatFlow
  {
  public:
	  /*! CBaseHeatFlow constructor
	  *
	  *	\param pMesh input closed mesh
	  * \param pWMesh sliced open mesh
	  */

    CBaseHeatFlow( M * pMesh, M * pWMesh );
	/*! CBaseHeatFlow destructor
	*/

    ~CBaseHeatFlow();

	/*!	Compute harmonic closed forms on the input mesh
	*/
    void calculate_harmonic_form();

  protected:
    /*! Compute the angle structure */
    virtual void _angle_structure() = 0;
	/*! Input closed mesh */	
    M * m_pMesh;
	/*! inpute mesh, sliced open along the shortest paths connecting two boundary loops */    
	M * m_pWMesh;

	typename M::CBoundary m_boundary;

	/*! Compute the harmonic 1-form */
    void _harmonic_1_form();
	/*! Integrate the closed form on the mesh */
    void _integrate();

  };

template<typename M>
CBaseHeatFlow<M>::CBaseHeatFlow( M * pMesh, M * pWMesh ):m_boundary( pWMesh )
{
  m_pMesh  = pMesh;
  m_pWMesh = pWMesh;
};

template<typename M>
CBaseHeatFlow<M>::~CBaseHeatFlow()
{
};


template<typename M>
void CBaseHeatFlow<M>::_harmonic_1_form()
{
    
    int vid = 0;

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
    {
		M::CVertex * v = *viter;
        v->idx() = vid ++;
    }

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
   {
	   M::CVertex * v = *viter;
      int  id = v->idx();
		
      double sum_b  = 0;

      for( M::VertexVertexIterator vviter( v ); !vviter.end();  ++vviter  )
      {
		  M::CVertex * w = *vviter;
		  M::CEdge   * e = m_pMesh->vertexEdge( v, w );

          if( m_pMesh->edgeVertex1( e ) != v )
          {
			  sum_b += -( e->weight() * e->du() );
          }
          else
          {
            sum_b +=  e->weight()* e->du() ;
          }
      }

	  v->u() = sum_b;
	 // if( fabs( sum_b)  > 0.25 )
     // printf("(%d %f) ",  id, sum_b );
  }

	CPoisson<M> P(m_pMesh);
	P.solve();

    for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
    {
		M::CEdge * e = *eiter;
		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );
		e->du() += v2->u() - v1->u();
	}


	 M::CBoundary boundary( m_pMesh );
 
	for( size_t i = 0; i < boundary.loops().size(); i ++ )
	{
		M::CLoop * pL = boundary.loops()[i];

		std::list<M::CHalfEdge*> pHs = pL->halfedges();
		double sum = 0; 
		for( std::list<M::CHalfEdge*>::iterator hiter = pHs.begin(); hiter != pHs.end(); hiter ++)
		  {
			  M::CHalfEdge * h = *hiter;
			  M::CEdge * e = m_pMesh->halfedgeEdge( h );
			  sum += (m_pMesh->edgeHalfedge(e,0) == h )? e->du():-e->du();
		  }
		  printf("Loop Integration %f\n", sum );
	}

}

template<typename M>
void CBaseHeatFlow<M>::_integrate()
{
	for( M::MeshEdgeIterator eiter( m_pWMesh ); !eiter.end(); eiter ++ )
  {
	  M::CEdge * we = *eiter;

	  M::CVertex * w1 = m_pWMesh->edgeVertex1( we );
	  M::CVertex * w2 = m_pWMesh->edgeVertex2( we );

	  //for diffuse exact form of multiply connected domain
	  //we use vertex id to replace vertex father

	  int id1 = ( w1->father() )? w1->father():w1->id();
	  int id2 = ( w2->father() )? w2->father():w2->id();

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
}

template<typename M>
void CBaseHeatFlow<M>::calculate_harmonic_form()
{
	_angle_structure();
	_harmonic_1_form();
    _integrate();
}

}	//namespace Holomorphy
}   //namespace MeshLib
#endif _BASE_HEAT_FLOW_