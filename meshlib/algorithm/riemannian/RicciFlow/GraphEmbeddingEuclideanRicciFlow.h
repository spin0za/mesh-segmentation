/*! \file GraphEmbeddingEuclideanRicciFlow.h
 *  \brief Euclidean Ricci flow algorithm for planar graph embedding
 *  \author David Gu
 *  \date   documented on 01/10/2014
 *
 *	Algorithm for Euclidean Ricci Flow for graph embedding
 */

#ifndef _GRAPH_EMBEDDING_EUCLIDEAN_RICCI_FLOW_H_
#define _GRAPH_EMBEDDING_EUCLIDEAN_RICCI_FLOW_H_

#include <map>
#include <vector>
#include <Eigen/Sparse>

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


namespace MeshLib
{

namespace RicciFlow
{

/*! \brief Class CGraphEmbeddingEuclideanRicciFlow
*
*	Algorithm for computing Euclidean Ricci flow graph embedding
*/
template<typename M>
class CGraphEmbeddingEuclideanRicciFlow : public CBaseRicciFlow<M>
  {
  public:
    /*! \brief CGraphEmbeddingEuclideanRicciFlow constructor
	 *  \param pMesh the input mesh
	 */
	  CGraphEmbeddingEuclideanRicciFlow( M * pMesh );
    /*! \brief CTangentialRicciFlow destructor
	 */
	  ~CGraphEmbeddingEuclideanRicciFlow(){};
	/*!	Computing the metric
	 */
	void _calculate_metric( double error_threshold, double step_length );

	/*!
	 *	Curvature flow, override 
	 */
     bool _flow( double );

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
	  *	Newton's method to optimize the entropy energy
	  * \param threshold err bound
	  * \param step_length step length
	  */
     virtual void   _Newton( double threshold, double step_length );

  };


template<typename M>
CGraphEmbeddingEuclideanRicciFlow<M>::CGraphEmbeddingEuclideanRicciFlow( M * pMesh ):CBaseRicciFlow( pMesh)
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		switch( pV->type() )
		{
		case 0: //vertex-vertex
			pV->u() = 0;
			break;
		case 1: //edge-vertex
			pV->u() = -1e+10;
			break;
		case 2: //face-vertex
			pV->u() = 0;
			break;
		};
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * pE = * eiter;
		M::CVertex * pV1 = m_pMesh->edgeVertex1( pE );
		M::CVertex * pV2 = m_pMesh->edgeVertex2( pE );

		if( ( pV1->type() == 0 && pV2->type() == 2 )|| ( pV1->type() == 2 && pV2->type() == 0 ) )
		{
			pE->inversive_distance() = 0;
		}
		else
		{
			pE->inversive_distance() = 1;
		}
	}

};

//Compute the edge length
template<typename M>
void CGraphEmbeddingEuclideanRicciFlow<M>::_length( double u1, double u2, typename M::CEdge * e )
{
	M::CVertex * pV1 = m_pMesh->edgeVertex1( e );
	M::CVertex * pV2 = m_pMesh->edgeVertex2( e );

	double r1 = exp( pV1->u() );
	double r2 = exp( pV2->u() );

	if( pV1->type() == 1 )
	{
		e->length() = r2;
		return;
	}

	if( pV2->type() == 1 )
	{
		e->length() = r1;
		return;
	}

	e->length() = sqrt( r1*r1 + r2*r2 );
	
};


//Calculate corner angle
template<typename M>
double CGraphEmbeddingEuclideanRicciFlow<M>::_cosine_law( double a, double b, double c )
{
          double cs =  ( a*a + b * b  - c * c )/( 2.0 * a * b );
          assert( cs <= 1.0 && cs >= -1.0 );
          return acos( cs );    
};


//Calculate edge weight

template<typename M>
void CGraphEmbeddingEuclideanRicciFlow<M>::_calculate_edge_weight()
{
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
  {
	  M::CEdge * e = *eiter;
      e->weight() = 0.0;
  }

	for(  M::MeshFaceIterator fiter( m_pMesh ) ; !fiter.end(); fiter ++ )
   {
	  M::CFace * f = *fiter;
	  M::CHalfEdge * head = NULL;

	  for( M::FaceHalfedgeIterator hiter(  f ); !hiter.end(); ++hiter )
      {
		  M::CHalfEdge * he = *hiter;
		  M::CVertex   * pS = m_pMesh->halfedgeSource( he );
		  M::CVertex   * pT = m_pMesh->halfedgeTarget( he );
		  if( pS->type() == 0 && pT->type() == 2 ) 
		  {
			  head = he;
			  break;
		  }
		  if( pS->type() == 2 && pT->type() == 0 ) 
		  {
			  head = he;
			  break;
		  }
      }
	  assert( head != NULL );
	  	  
	  M::CHalfEdge * ph = m_pMesh->halfedgePrev( head );
	  M::CHalfEdge * nh = m_pMesh->halfedgeNext( head );

	  M::CEdge * pe = m_pMesh->halfedgeEdge( ph   );
	  M::CEdge * ne = m_pMesh->halfedgeEdge( nh   );
	  M::CEdge * pE = m_pMesh->halfedgeEdge( head );

	  pe->weight() = 0;
	  ne->weight() = 0;
	
	  double h = pe->length() * ne->length() / pE->length();

	  pE->weight() += h/pE->length();

  }
};

//set target curvature

template<typename M>
void CGraphEmbeddingEuclideanRicciFlow<M>::_set_target_curvature()
{
};

//compute metric

template<typename M>
void CGraphEmbeddingEuclideanRicciFlow<M>::_calculate_metric( double error_threshold, double step_length )
{
	_Newton( error_threshold, step_length );
};

//gradient flow method

template<typename M>
bool CGraphEmbeddingEuclideanRicciFlow<M>::_flow( double error_threshold )
{
  int num = m_pMesh->numVertices();

  for( int k = 0; k < 64000; k ++  )
	  {
 		  _calculate_edge_length();
	      _set_target_curvature();
		  _calculate_edge_weight();

		  _calculate_corner_angle();
		  _calculate_vertex_curvature();

		  double error =  _calculate_curvature_error();
		  printf("Current error is %f\r\n", error );

		  if( error < error_threshold)  return true;
  	  

	  //set b 
		  for(M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
		  {
			  M::CVertex * v = *viter;
			  if( v->type() == 1 ) continue; //edge vertex
			  double dif = v->target_k() - v->k();
			  v->u() += dif * 5e-2;
		  }
    }
    return false;
};


//Newton's method for optimizing entropy energy

template<typename M>
void CGraphEmbeddingEuclideanRicciFlow<M>::_Newton( double threshold, double step_length )
{
	int num = m_pMesh->numVertices();

	
  	while( true )
	{
		//the order of the following functions really matters

 		_calculate_edge_length();
		_calculate_corner_angle();
		_calculate_vertex_curvature();
		_calculate_edge_weight();

		double error =  _calculate_curvature_error();
		printf("Newton's Method: Current error is %f\r\n", error );
		if( error < threshold) break;
	

		for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
		{
			M::CVertex * v = *viter;
			if( v->type() == 1 ) continue;
			v->huv()[0] = v->u();	
			v->u() = v->target_k() - v->k();
		}

		CPoisson<M> P(m_pMesh );
		P.solve3();
	

		for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
		{
			M::CVertex * v = *viter;
			if( v->type() == 1 ) continue;
			v->u() = v->huv()[0] + v->u() * step_length;
		}

  }

};


} //namespace RicciFlow

} //namespace MeshLib	

#endif  _GRAPH_EMBEDDING_EUCLIDEAN_RICCI_FLOW_H_