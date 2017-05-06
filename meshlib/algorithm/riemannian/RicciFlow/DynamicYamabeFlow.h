/*! \file   DynamicYamabeFlow.h
 *  \brief  Dynamic Euclidean Yamabe flow
 *  \author David Gu
 *  \date   documented on 02/08/2013
 *
 *	Algorithm for Dynamic Yamabe Flow
 */

#ifndef _DYNAMIC_YAMABE_FLOW_H_
#define _DYNAMIC_YAMABE_FLOW_H_

#include <map>
#include <vector>

#include "Geometry/Circle.h"
#include "Mesh/DynamicMesh.h"
#include "DynamicRicciFlowMesh.h"
#include "DynamicBaseRicciFlow.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "mesh/boundary.h"
#include "Parser/parser.h"
#include "Parser/traits_io.h"
#include "Geometry/Circle.h"

namespace MeshLib
{

namespace RicciFlow
{
	/*! \brief CDynamicYamabeFlow class
	 *  
	 *  Algorithm for Dynamic Euclidean Yamabe flow
	 */
	 
 template<typename M>
 class CDynamicYamabeFlow : public CDynamicBaseRicciFlow<M>
  {
  public:
	  /*! \brief CDynamicYamabeFlow constructor
	   *
	   *  call base class constructor 
	   */
	  CDynamicYamabeFlow( M * pMesh );
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
	  *	Normalization
	  * \param du the du vector
	  * \param n dimension of the du vector
	  */
	void _normalization( Eigen::VectorXd & du, int n );

	/*!
	 *	Calculate the edge weight
	 */
    void _calculate_edge_weight();

	/*!
	 *	Set the target curvature on each vertex
	 */
    virtual void    _set_target_curvature();

	/*!
	 *	Set initial edge length
	 */
	virtual void _calculate_initial_edge_length();
	/*!
	 *	Preserve Delaunay	
	 */
	void _preserve_Delaunay();
  
 };

template<typename M>
CDynamicYamabeFlow<M>::CDynamicYamabeFlow( M * pMesh ): CDynamicBaseRicciFlow( pMesh)
{
	_calculate_initial_edge_length();
}

 //calculate edge length
template<typename M>
void CDynamicYamabeFlow<M>::_length( double u1, double u2, typename M::CEdge * e )
{
	  e->length() = exp(u1 ) * exp( u2 ) * e->initial_length();
};


 //calculate edge weight

template<typename M>
void CDynamicYamabeFlow<M>::_calculate_initial_edge_length()
{

	  for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	  {
		  M::CEdge * e = *eiter;
		  e->initial_length() = m_pMesh->edgeLength( e );
	  }
}


 //calculate edge weight

template<typename M>
void CDynamicYamabeFlow<M>::_calculate_edge_weight()
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
void CDynamicYamabeFlow<M>::_calculate_metric()
{

	  _calculate_initial_edge_length();

	  _set_target_curvature();
	  double error_threshold = 1e-8;
	  //double step_length = 0.5;
	  double step_length = 1.0;
	  //_Newton( error_threshold, step_length );
	  //_flow( error_threshold );	  
	  _Dynamic_Newton( error_threshold, step_length );
};      


//Euclidean cosine law
template<typename M>
double CDynamicYamabeFlow<M>::_cosine_law( double a, double b, double c )
{
          double cs =  ( a*a + b * b  - c * c )/( 2.0 * a * b );
          assert( cs <= 1.0 && cs >= -1.0 );
          return acos( cs );    
};

//normalization

template<typename M>
void CDynamicYamabeFlow<M>::_normalization( Eigen::VectorXd & x, int num )
{

	double s = 0;
	for(int i = 0; i < num; i ++ )
	{
		s += x(i);
	}
	s /= num;

	for (int i = 0; i < num; i++)
	{
	 x(i) -= s;
	}
};


//set target curvature

template<typename M>
void CDynamicYamabeFlow<M>::_set_target_curvature()
{
	int Euler = m_pMesh->numVertices() + m_pMesh->numFaces() - m_pMesh->numEdges();
	double total_curvature = 2 * PI * Euler;
/*
	//extremal curvature
	double w = 0;
	for( CDynamiM::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		V * v = *viter;
		v->target_k() = 2*PI-0.1;
		w += v->target_k();
	}
	
  for( CDynamiM::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
    V * v = *viter;
    v->target_k() += total_curvature - w;
	break;
  }


	//uniformization metric
  for( CDynamiM::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
    V * v = *viter;
    v->target_k() = total_curvature/m_pMesh->numVertices();
  }


	//uniformization metric
  for( CDynamiM::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
    V * v = *viter;
    v->target_k() = 0;
  }


  _read_vertex_target_k<CDynamiM, V,E,F,H>( m_pMesh );
*/

	//uniformization metric
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
	  v->target_k() = 0;
  }

};



//preserve Delaunay Triangulation

template<typename M>
void CDynamicYamabeFlow<M>::_preserve_Delaunay()
{
	m_pMesh->_preserve_Delaunay();

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );
		e->initial_length() = e->length()/exp(v1->u()+v2->u());
	}
};





} //namespace RicciFlow
} //namespace MeshLib

#endif  _DYNAMIC_YAMABE_FLOW_H_