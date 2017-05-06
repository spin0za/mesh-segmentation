/*! \file HyperbolicYamabeFlow.h
 *  \brief Hyperbolic Yamabe flow algorithm
 *  \author David Gu
 *  \date   documented on 11/14/2010
 *
 *	Algorithm for Inversive Distance Hyperbolic Ricci Flow
 */

#ifndef _HYPERBOLIC_YAMABE_FLOW_H_
#define _HYPERBOLIC_YAMABE_FLOW_H_

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
#include "Geometry/Hyperbolic.h"

namespace MeshLib
{

namespace RicciFlow
{
	
/*! \brief Class CHyperbolicYamabeFlow
*
*	Algorithm for computing Hyperbolic Yamabe flow
*/
template<typename M>
class CHyperbolicYamabeFlow : public CBaseHyperbolicRicciFlow<M>
  {
  public:
    /*! \brief CHyperbolicYamabeFlow constructor
	 *  \param pMesh the input mesh
	 */
	  CHyperbolicYamabeFlow( M * pMesh );
    /*! \brief CHyperbolicYamabeFlow destructor
	 */
	  ~CHyperbolicYamabeFlow(){};

  protected:

		/*!	Computing the edge weight
		 */
		void _calculate_edge_weight();
	

	  /*!
	   *	Calculate each edge length, has to be defined in the derivated classes
	   */
	  void _length( double u1, double u2, typename M::CEdge * e );

	/*!
	 *	Set the target curvature on each vertex
	 */
    void    _set_target_curvature();


  private:
	   /*!
	    *	Compute the \f$\frac{\partial \theta_i }{\partial u_j}\f$
		*/
	   void _compute_half_edge_weight_ij();
	   /*!
	    *	Compute the \f$\frac{\partial \theta_i }{\partial u_i}\f$
		*/
	   void _compute_half_edge_weight_ii();
	   /*!
	    *	Compute the \f$\frac{\partial \theta_i }{\partial u_j}\f$
		*   \param y_i,y_j,y_k edge lengths
		*/

	   double _dtheta_i_du_j( double y_i, double y_j, double y_k );

	   /*!
	    *	Compute the \f$\frac{\partial \theta_i }{\partial u_i}\f$
		*   \param y_i,y_j,y_k edge lengths
		*/
	   double _dtheta_i_du_i( double y_i, double y_j, double y_k );

  };




template<typename M>
double CHyperbolicYamabeFlow<M>::_dtheta_i_du_i( double y_i, double y_j, double y_k )
{
	double c_i = cosh( y_i );
	double c_j = cosh( y_j );
	double c_k = cosh( y_k );

	double theta_i = _cosine_law( y_j, y_k, y_i );
	double A = sin(theta_i) * sinh(y_j) * sinh(y_k);
	A = 1/A;

	double d = -A * ( 2 * c_i * c_j * c_k - c_j*c_j - c_k * c_k + c_i * c_j + c_i * c_k - c_j - c_k )/( (c_j + 1) *( c_k +1) );

	return d;
};

template<typename M>
double CHyperbolicYamabeFlow<M>::_dtheta_i_du_j( double y_i, double y_j, double y_k )
{
	double c_i = cosh( y_i );
	double c_j = cosh( y_j );
	double c_k = cosh( y_k );

	double theta_i = _cosine_law(  y_j, y_k, y_i );
	double A = sin(theta_i) * sinh(y_j) * sinh(y_k);
	A = 1/A;

	double d = A * ( c_i + c_j - c_k - 1 )/( c_k  + 1 );
	return d;
};

template<typename M>
void CHyperbolicYamabeFlow<M>::_compute_half_edge_weight_ij( )
{
	for(  M::MeshHalfEdgeIterator hiter( m_pMesh ); !hiter.end(); hiter ++ )
	{
		M::CHalfEdge * pH = *hiter;
		M::CHalfEdge * pN = m_pMesh->faceNextCcwHalfEdge( pH );
		M::CHalfEdge * pR = m_pMesh->faceNextClwHalfEdge( pH );

		double y_i = m_pMesh->halfedgeEdge( pR )->length();
		double y_j = m_pMesh->halfedgeEdge( pN )->length();
		double y_k = m_pMesh->halfedgeEdge( pH )->length();

		pH->c_w_ij() = -_dtheta_i_du_j( y_i, y_j, y_k );
	}
};

template<typename M>
void CHyperbolicYamabeFlow<M>::_compute_half_edge_weight_ii()
{
	for( M::MeshHalfEdgeIterator hiter( m_pMesh ); !hiter.end(); hiter ++ )
	{
		M::CHalfEdge * pH = *hiter;
		M::CHalfEdge * pN = m_pMesh->faceNextCcwHalfEdge( pH );
		M::CHalfEdge * pR = m_pMesh->faceNextClwHalfEdge( pH );

		double y_i = m_pMesh->halfedgeEdge( pR )->length();
		double y_j = m_pMesh->halfedgeEdge( pN )->length();
		double y_k = m_pMesh->halfedgeEdge( pH )->length();

		pH->c_w_ii() = -_dtheta_i_du_i( y_i, y_j, y_k );
	}
};



template<typename M>
CHyperbolicYamabeFlow<M>::CHyperbolicYamabeFlow( M * pMesh ):CBaseHyperbolicRicciFlow( pMesh)
{
	int vid = 0;
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		pV->idx() = vid ++;
	}
};

//Compute the edge length
template<typename M>
void CHyperbolicYamabeFlow<M>::_length( double u1, double u2, typename M::CEdge * e )
{
	//key
	double l = _asinh( sinh( m_pMesh->edgeLength( e )/2.0 ) * exp( u1 + u2 ) ) * 2.0;
    e->length() = l;
};

//set target curvature

template<typename M>
void CHyperbolicYamabeFlow<M>::_set_target_curvature()
{
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
   M::CVertex* v = *viter;
    v->target_k() = 0;
  }
};


//Calculate edge weight

template<typename M>
void CHyperbolicYamabeFlow<M>::_calculate_edge_weight()
{
	_compute_half_edge_weight_ij();
	_compute_half_edge_weight_ii();

};



} //namespace RicciFlow
} //namespace MeshLib
#endif  _HYPERBOLIC_YAMABE_FLOW_H_