/*! \file YamabeFlow.h
 * \brief Euclidean Yamabe flow
 *  \author David Gu
 *  \date   documented on 11/13/2010
 *
 *	Algorithm for Yamabe Flow
 */

#ifndef _YAMABE_FLOW_H_
#define _YAMABE_FLOW_H_

#include <map>
#include <vector>

#include "Mesh/BaseMesh.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "mesh/boundary.h"
#include "Parser/parser.h"
#include "RicciFlowMesh.h"
#include "BaseRicciFlow.h"
#include "Geometry/Circle.h"

namespace MeshLib
{

namespace RicciFlow
{
	/*! \brief CYamabeFlow class
	 *  
	 *  Algorithm for Euclidean Yamabe flow
	 */
	 
 template<typename M>
 class CYamabeFlow : public CBaseRicciFlow<M>
  {
  public:
	  /*! \brief CYamabeFlow constructor
	   *
	   *  call base class constructor 
	   */
	  CYamabeFlow( M * pMesh );
	  /*! override the virtual function calculate_metric */
	  virtual void _calculate_metric();
	

  protected:

	  /*!
	   *	Calculate each edge length, has to be defined in the derivated classes
	   */
	  void _length( double u1, double u2, typename M::CEdge * e );

	/*!
	 *	Cosine law, has to be defined in the derivated classes
	 */
	double _cosine_law( double a, double b, double c ) ;


	/*!
	 *	Calculate the edge weight
	 */
    void _calculate_edge_weight();

	/*!
	 *	Set the target curvature on each vertex
	 */
    virtual void    _set_target_curvature();

  };

template<typename M>
CYamabeFlow<M>::CYamabeFlow( M * pMesh ): CBaseRicciFlow( pMesh)
{
}

 //calculate edge length
template<typename M>
void CYamabeFlow<M>::_length( double u1, double u2, typename M::CEdge * e )
{
	  double l  = m_pMesh->edgeLength( e );
	  e->length() = exp(u1 ) * exp( u2 ) * l;
};

 //calculate edge weight

template<typename M>
void CYamabeFlow<M>::_calculate_edge_weight()
{

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	  {
		  M::CEdge * e = *eiter;
		  e->weight() = 0.0;

		  M::CHalfEdge * h = m_pMesh->edgeHalfedge( e, 0 );
		  M::CHalfEdge * n = m_pMesh->faceNextCcwHalfEdge( h );
		  e->weight() += cos( n->angle() )/sin( n->angle() );

		  h = m_pMesh->edgeHalfedge( e, 1 );
		  if( h == NULL ) continue;
		  n = m_pMesh->faceNextCcwHalfEdge( h );
		  e->weight() += cos( n->angle() )/sin( n->angle() );
	  }
}

template<typename M>
void CYamabeFlow<M>::_calculate_metric()
{

	  _set_target_curvature();
	  double error_threshold = 1e-6;
	  //double step_length = 0.5;
	  double step_length = 0.5;
	  _Newton( error_threshold, step_length );
		  
};      


//Euclidean cosine law
template<typename M>
double CYamabeFlow<M>::_cosine_law( double a, double b, double c )
{
          double cs =  ( a*a + b * b  - c * c )/( 2.0 * a * b );
          assert( cs <= 1.0 && cs >= -1.0 );
          return acos( cs );    
};


//set target curvature

template<typename M>
void CYamabeFlow<M>::_set_target_curvature()
{
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
	  v->target_k() = 0;
  }
};

}	//namespace RicciFlow
}	//namespace MeshLib

#endif  _YAMABE_FLOW_H_