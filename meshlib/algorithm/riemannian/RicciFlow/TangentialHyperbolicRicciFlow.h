/*! \file TangentialHyperbolicRicciFlow.h
 *  \brief Tangential Hyperbolic Ricci flow algorithm
 *  \author David Gu
 *  \date   documented on 10/17/2010
 *
 *	Algorithm for general Ricci Flow
 */

#ifndef _TANGENTIAL_HYPERBOLIC_RICCI_FLOW_H_
#define _TANGENTIAL_HYPERBOLIC_RICCI_FLOW_H_

#include <map>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "Mesh/BaseMesh.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "mesh/boundary.h"
#include "Parser/parser.h"
#include "RicciFlowMesh.h"
#include "BaseHyperbolicRicciFlow.h"


namespace MeshLib
{
namespace RicciFlow
{
	
	
/*! \brief Class CTangentialHyperbolicRicciFlow
*
*	Algorithm for computing Hyperbolic Ricci flow using tangential Circle packing metric
*/
template<typename M>
class CTangentialHyperbolicRicciFlow : public CBaseHyperbolicRicciFlow<M>
  {
  public:
    /*! \brief CTangentialRicciFlow constructor
	 *  \param pMesh the input mesh
	 */
	  CTangentialHyperbolicRicciFlow( M * pMesh );
    /*! \brief CTangentialRicciFlow destructor
	 */
	  ~CTangentialHyperbolicRicciFlow(){};
	

  protected:

	  /*!
	   *	Calculate each edge length, has to be defined in the derivated classes
	   */
	  void _length( double u1, double u2, typename M::CEdge * e );


	/*!
	 *	Calculate the edge weight
	 */

	void _calculate_edge_weight();
	
	/*!
	 *	Calculate the edge weight for a halfedge
	 *  \param pH input halfedge
	 */
	void _calculate_edge_weight( typename M::CHalfEdge * pH );

	/*!
	 *	Set the target curvature on each vertex
	 */
    virtual void    _set_target_curvature();
	
  
  };

template<typename M>
void CTangentialHyperbolicRicciFlow<M>::_calculate_edge_weight()
{
  for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); ++ fiter)
  {
	  M::CFace * f = *fiter;
	  for( M::FaceHalfedgeIterator hiter( f ); !hiter.end(); hiter ++ )
	  {
		  M::CHalfEdge * h = *hiter;
	      _calculate_edge_weight( h );
	  }
  }
};


template<typename M>
CTangentialHyperbolicRicciFlow<M>::CTangentialHyperbolicRicciFlow( M * pMesh ) : CBaseHyperbolicRicciFlow<M>( pMesh )
{

	//initialize u value at each vertex
	
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * pV = *viter;
		pV->u() = _r2u( 1.0 );
	}
};

//Compute the edge length
template<typename M>
void CTangentialHyperbolicRicciFlow<M>::_length( double u1, double u2, typename M::CEdge * e )
{
	double r1 = _u2r( u1 );
    double r2 = _u2r( u2 );
    double l = _inverse_cosine_law( r1,r2, PI );
	e->length() = l;
	e->length() = r1 + r2;
};

//set target curvature

template<typename M>
void CTangentialHyperbolicRicciFlow<M>::_set_target_curvature()
{
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
	  v->target_k() = 0;
  }
};


template<typename M>
void CTangentialHyperbolicRicciFlow<M>::_calculate_edge_weight( typename M::CHalfEdge * pH )
{
	M::CHalfEdge * e_k = pH;
    M::CHalfEdge * e_i = m_pMesh->faceNextCcwHalfEdge( e_k );
    M::CHalfEdge * e_j = m_pMesh->faceNextCcwHalfEdge( e_i );

	M::CVertex * v_i = m_pMesh->halfedgeSource( e_k );
	M::CVertex * v_j = m_pMesh->halfedgeSource( e_i );
	M::CVertex * v_k = m_pMesh->halfedgeSource( e_j );


    double r_i = _u2r( v_i->u() );
    double r_j = _u2r( v_j->u() );
    double r_k = _u2r( v_k->u() );

    double theta_i = e_j->angle();
    double theta_j = e_k->angle();
    double theta_k = e_i->angle();

    double x_i = m_pMesh->halfedgeEdge( e_i )->length();
    double x_j = m_pMesh->halfedgeEdge( e_j )->length();
    double x_k = m_pMesh->halfedgeEdge( e_k )->length();

    double A_ijk = sinh(x_i) * sinh(x_j) * sin( theta_k); 

    /*!
	* \f$ d_ij = \frac{\partial \theta_i}{\partial x_j} \f$
	*/	
    double d_ii =   sinh( x_i ) / A_ijk;
    double d_ij = - d_ii * cos( theta_k );
    double d_ik = - d_ii * cos( theta_j );

    double w_ii = d_ij + d_ik;
    double w_ij = d_ii + d_ik;
    double w_ik = d_ii + d_ij;
    
    /*!
	 * \f$ w_ij = \frac{\partial \theta_i}{\partial r_j} * s(r_j) \f$
	 */
    double D_ij = sinh( r_j ) * w_ij;
    double D_ii = sinh( r_i ) * w_ii;
    double D_ik = sinh( r_k ) * w_ik;

    //e_k->c_w_ii() = -D_ii;
    //e_k->c_w_ij() = -D_ij;

	M::CHalfEdge * sh_k = m_pMesh->halfedgeSym( e_k );

    sh_k->c_w_ii() = -D_ii;
    sh_k->c_w_ij() = -D_ij;
	
};

}
}

#endif  _TANGENTIAL_HYPERBOLIC_RICCI_FLOW_H_