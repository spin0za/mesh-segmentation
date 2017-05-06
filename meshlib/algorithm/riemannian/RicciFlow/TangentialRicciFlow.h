/*! \file TangentialRicciFlow.h
 *  \brief General Euclidean Ricci flow algorithm
 *  \author David Gu
 *  \date   documented on 10/17/2010
 *
 *	Algorithm for general Ricci Flow
 */

#ifndef _TANGENTIAL_RICCI_FLOW_H_
#define _TANGENTIAL_RICCI_FLOW_H_

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

/*! \brief Class CTangentialRicciFlow
*
*	Algorithm for computing Ricci flow using Tangential Circle Packing
*/
template<typename M>
class CTangentialRicciFlow : public CBaseRicciFlow<M>
  {
  public:
    /*! \brief CTangentialRicciFlow constructor
	 *  \param pMesh the input mesh
	 */
	  CTangentialRicciFlow( M * pMesh );
    /*! \brief CTangentialRicciFlow destructor
	 */
	  ~CTangentialRicciFlow(){};
	/*!	Computing the metric
	 */
	void _calculate_metric();
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

  };


template<typename M>
CTangentialRicciFlow<M>::CTangentialRicciFlow( M * pMesh ):CBaseRicciFlow( pMesh)
{
};

//Compute the edge length
template<typename M>
void CTangentialRicciFlow<M>::_length( double u1, double u2, typename M::CEdge * e )
{
	  e->length() = exp(u1) + exp(u2);
};


//Calculate corner angle
template<typename M>
double CTangentialRicciFlow<M>::_cosine_law( double a, double b, double c )
{
          double cs =  ( a*a + b * b  - c * c )/( 2.0 * a * b );
          assert( cs <= 1.0 && cs >= -1.0 );
          return acos( cs );    
};


//Calculate edge weight

template<typename M>
void CTangentialRicciFlow<M>::_calculate_edge_weight()
{
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
  {
	  M::CEdge * e = *eiter;
      e->weight() = 0.0;
  }

	for(  M::MeshFaceIterator fiter( m_pMesh ) ; !fiter.end(); fiter ++ )
  {
	  M::CFace * f = *fiter;

      double r[3];
      int i = 0;
	  for( M::FaceHalfedgeIterator hiter(  f ); !hiter.end(); ++hiter )
      {
		  M::CHalfEdge * he = *hiter;
		  M::CVertex * v = m_pMesh->halfedgeTarget(he);
         r[i++] = exp( v->u() );
      }

      double w = sqrt(r[0]*r[1]*r[2]/(r[0]+r[1]+r[2]));
      
	  for( M::FaceEdgeIterator eiter(f); !eiter.end(); ++ eiter )
      {
		  M::CEdge * e = * eiter;
          e->weight() += w/e->length();
      }
  }
};

//set target curvature

template<typename M>
void CTangentialRicciFlow<M>::_set_target_curvature()
{
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
    v->target_k() = 0;
  }

  std::vector<M::CLoop*>& pLs = m_boundary.loops();
	
  int id = 0;

  for( std::vector<M::CLoop*>::iterator liter = pLs.begin(); liter != pLs.end(); liter++)
  {
	  M::CLoop * pL = *liter;
	  std::list<M::CHalfEdge*> & pHes = pL->halfedges();

	  double sum = 0;
	  double inv_sum = 0;
	
	  for( std::list<M::CHalfEdge*>::iterator hiter = pHes.begin(); hiter != pHes.end(); hiter ++ )
	  {
		  M::CHalfEdge  * he = *hiter;
		  M::CEdge  * pe = m_pMesh->halfedgeEdge( he );

			sum  += pe->length();
			inv_sum += 1/pe->length();
 	  }

	  for( std::list<M::CHalfEdge*>::iterator hiter = pHes.begin(); hiter != pHes.end(); hiter ++ )
	  {
		  M::CHalfEdge * ce = *hiter;
		  M::CVertex   * pv = m_pMesh->halfedgeTarget( ce );
		  M::CHalfEdge * he = m_pMesh->vertexMostCcwInHalfEdge( pv );
		  M::CHalfEdge * te =  m_pMesh->vertexMostClwOutHalfEdge( pv );

			double L = ( m_pMesh->halfedgeEdge(he)->length() + m_pMesh->halfedgeEdge( te )->length()  )/2.0;

		 // map all boundaries to circular holes
			if( id == 0 )
			{
				double tk = 2*PI*L/sum;
				pv->target_k() = tk;
			}
			else
			{
				double tk = -2*PI*L/sum;
				pv->target_k() = tk;
			}
	  }
	  id ++;
  }

};

//compute metric

template<typename M>
void CTangentialRicciFlow<M>::_calculate_metric()
{

  //double error = 1e-6;
  double error = 5e-4;

  _calculate_edge_length();

  while( true )
  {
	  _set_target_curvature();
	  _Newton( error, 1 );
    //break;
      if( _flow( error ) ) break;
}      


};

//gradient flow method

template<typename M>
bool CTangentialRicciFlow<M>::_flow( double error_threshold )
{
  int num = m_pMesh->numVertices();

  for( int k = 0; k < 64; k ++  )
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
			double dif = v->target_k() - v->k();
            v->u() += dif * 2e-2;
		  }
    }
    return false;
};

} //namespace RicciFlow

} //namespace MeshLib	

#endif  _TANGENTIAL_RICCI_FLOW_H_