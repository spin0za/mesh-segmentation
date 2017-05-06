/*! \file TangentialMixedHyperbolicRicciFlow.h
 *  \brief Tangential Mixed Hyperbolic Ricci flow algorithm
 *  \author David Gu
 *  \date   documented on 08/03/2011
 *
 *	Algorithm for tangential mixed hyperbolic Ricci Flow
 */

#ifndef _TANGENTIAL_MIXED_HYPERBOLIC_RICCI_FLOW_H_
#define _TANGENTIAL_MIXED_HYPERBOLIC_RICCI_FLOW_H_
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
//#include "BaseRicciFlow.h"
#include "BaseMixedRicciFlow.h"
#include "TangentialHyperbolicRicciFlow.h"

namespace MeshLib
{
namespace RicciFlow
{
	
/*! \brief Class CTangentialHyperbolicRicciFlow
*
*	Algorithm for computing Hyperbolic Ricci flow using tangential Circle packing metric with mixed target curvature on interior vertices and target conformal factor on boundary vertices
*/
template<typename M>
class CTangentialMixedHyperbolicRicciFlow : public CBaseMixedRicciFlow<M>
//class CTangentialMixedHyperbolicRicciFlow : public CBaseRicciFlow<M>
  {
  public:
    /*! \brief CTangentialMixedRicciFlow constructor
	 *  \param pMesh the input mesh
	 */
	  CTangentialMixedHyperbolicRicciFlow( M * pMesh );
    /*! \brief CTangentialMixedRicciFlow destructor
	 */
	  ~CTangentialMixedHyperbolicRicciFlow(){};
	/*!	Computing the metric
	 */
	void _calculate_metric();
	

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
	 *	Cosine law, has to be defined in the derivated classes
	 */
	double _inverse_cosine_law( double a, double b, double c ) ;
	 /*!
	  *	Normalization
	  * \param du the du vector
	  * \param n dimension of the du vector
	  * do nothing
	  */
	void _normalization( Eigen::VectorXd & du, int n );


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
	 *	Set the target curvature on each interior vertex
	 */
    virtual void    _set_target_curvature();
	/*!
	 *	Set the target conformal factor on each boundary vertex
	 */
    virtual void    _set_target_conformal_factor();
	/*!
	 *	Curvature flow 
	 */
    bool   _flow( double );

	 /*!
	  *	calculate hessian matrix Hessain 
	  * \param SparseMatrix
	  */
	void _calculate_Hessain( Eigen::SparseMatrix<double> & M );

  protected:
	  /*!
	   *	converting u to radius r
	   */
	  double _u2r( double u );
	  /*! 
	   *	converting radius r to u
	   */
	  double _r2u( double r );
  };


template<typename M>
CTangentialMixedHyperbolicRicciFlow<M>::CTangentialMixedHyperbolicRicciFlow( M * pMesh ):CBaseMixedRicciFlow( pMesh )
{
	//initialize u value at each vertex

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * pV = *viter;
		pV->u() = _r2u( 1.0 );
	}
};



//Calculate corner angle
template<typename M>
double CTangentialMixedHyperbolicRicciFlow<M>::_cosine_law( double a, double b, double c )
{
	double C;
	C = acos( (cosh(a) * cosh(b)-cosh(c))/(sinh(a)*sinh(b)) );
	return C;
};

//given edge lengths a,b and angle C, return edge length c

template<typename M>
double CTangentialMixedHyperbolicRicciFlow<M>::_inverse_cosine_law( double a, double b, double C )
{
  double c;
  c = cosh(a) * cosh(b) - sinh(a)*sinh(b) * cos(C) ;
  c = log( c+sqrt(c*c-1) );
  return c;
}

// u = log tanh (r/2)
template<typename M>
double CTangentialMixedHyperbolicRicciFlow<M>::_u2r( double u )
{  
    double e = exp(u);
    double r;
    r = log(  (1+e)/(1 - e ) );
    return r;
}

// u = log tanh (r/2)
template<typename M>
double CTangentialMixedHyperbolicRicciFlow<M>::_r2u( double r )
{  
    return log(tanh(r/2));
};

//calculate edge length
template<typename M>
void CTangentialMixedHyperbolicRicciFlow<M>::_length( double u1, double u2, typename M::CEdge * e )
{
	double radius_1 = _u2r( u1 );
	double radius_2 = _u2r( u2 );
	e->length() = radius_1 + radius_2;
}


//Calculate edge weight

template<typename M>
void CTangentialMixedHyperbolicRicciFlow<M>::_calculate_edge_weight()
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

//set target curvature

template<typename M>
void CTangentialMixedHyperbolicRicciFlow<M>::_set_target_curvature()
{
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
	if( m_pMesh->isBoundary( v ) ) continue;
    v->target_k() = 0;
  }
};

//set target conformal factor

template<typename M>
void CTangentialMixedHyperbolicRicciFlow<M>::_set_target_conformal_factor()
{
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
	if( !m_pMesh->isBoundary( v ) ) continue;
    v->u() = 0.0;
  }
};

//compute metric

template<typename M>
void CTangentialMixedHyperbolicRicciFlow<M>::_calculate_metric()
{

	_set_target_curvature();
	_set_target_conformal_factor();
	double error_threshold = 1e-6;
	double step_length = 0.1;
	_Newton( error_threshold, step_length );

};

//gradient flow method

template<typename M>
bool CTangentialMixedHyperbolicRicciFlow<M>::_flow( double error_threshold )
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

//normalization

template<typename M>
void CTangentialMixedHyperbolicRicciFlow<M>::_normalization( Eigen::VectorXd & x, int num )
{
	//for hyperbolic Ricci flow, do nothing
};


template<typename M>
void CTangentialMixedHyperbolicRicciFlow<M>::_calculate_edge_weight( typename M::CHalfEdge * pH )
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

    e_k->c_w_ii() = D_ii;
    e_k->c_w_ij() = D_ij;
};

template<typename M>
void CTangentialMixedHyperbolicRicciFlow<M>::_calculate_Hessain( Eigen::SparseMatrix<double> & M )
{
	std::vector<Eigen::Triplet<double> > _coefficients;

  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
  {
      M::CVertex * v = *viter;
	  if( m_pMesh->isBoundary( v ) ) continue;

      double a_ii = 0;

	  for( M::VertexOutHalfedgeIterator hiter( m_pMesh, v ); !hiter.end(); ++ hiter )
	  {
		  M::CHalfEdge * he = *hiter;
		M::CVertex * w = m_pMesh->halfedgeTarget( he );
		if( m_pMesh->isBoundary( w ) ) continue;
	  	a_ii += he->c_w_ii();
	  }
      //pMatrix->AddElementTail( v->idx(), v->idx(), -a_ii); 
	  _coefficients.push_back( Eigen::Triplet<double>( v->idx(), v->idx(), -a_ii) );
  }

  for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); ++ eiter )
  {
	  M::CEdge * e = *eiter;

      double a_ij = 0;

	  M::CHalfEdge * he = m_pMesh->edgeHalfedge(e,0);

      a_ij += he->c_w_ij();

	  M::CHalfEdge *  twin = m_pMesh->halfedgeSym( he );
      if( twin != NULL )
      {
            a_ij += twin->c_w_ij();
      }

      M::CVertex * v0 = m_pMesh->edgeVertex1( e );
      M::CVertex * v1 = m_pMesh->edgeVertex2( e );

	  if( m_pMesh->isBoundary( v0 ) || m_pMesh->isBoundary( v1 ) ) continue;

      //pMatrix->AddElementTail( v0->idx(), v1->idx(), -a_ij);   
      //pMatrix->AddElementTail( v1->idx(), v0->idx(), -a_ij);   
	  _coefficients.push_back( Eigen::Triplet<double>( v0->idx(), v1->idx(), -a_ij) );
	  _coefficients.push_back( Eigen::Triplet<double>( v1->idx(), v0->idx(), -a_ij) );
  }

  M.setFromTriplets( _coefficients.begin(), _coefficients.end() );
};

}
}

#endif  _TANGENTIAL_MIXED_HYPERBOLIC_RICCI_FLOW_H_