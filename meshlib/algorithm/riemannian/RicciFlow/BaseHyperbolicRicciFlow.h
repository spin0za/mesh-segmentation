/*! \file BaseHyperbolicRicciFlow.h
 *  \brief Base Hyperbolic Ricci flow algorithm
 *  \author David Gu
 *  \date   documented on 12/14/2013
 *
 *	Algorithm for general Hyperbolic Ricci Flow
 */

#ifndef _BASE_HYPERBOLIC_RICCI_FLOW_H_
#define _BASE_HYPERBOLIC_RICCI_FLOW_H_

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
#include "BaseRicciFlow.h"


namespace MeshLib
{
	
namespace RicciFlow
{

/*! 
*	\class CBaseHyperbolicRicciFlow
*	\brief Class CBaseHyperbolicRicciFlow
*
*	Algorithm for computing Hyperbolic Ricci flow
*/
template<typename M>
class CBaseHyperbolicRicciFlow : public CBaseRicciFlow<M>
  {
  public:
    /*! \brief CBaseRicciFlow constructor
	 *  \param pMesh the input mesh
	 */
	  CBaseHyperbolicRicciFlow( M * pMesh );
    /*! \brief CBaseRicciFlow destructor
	 */
	  ~CBaseHyperbolicRicciFlow(){};
	/*!	\brief Computing the metric
	 *  \param threshould  : error threshold
	 *  \param step_length : step length in the optimization
	 */
	void _calculate_metric(   double threshold = 1e-6, double step_length = 0.5 );

  protected:

	  /*!
	   *	\brief Calculate each edge length, has to be defined in the derivated classes
	   *    \param u1, u2 : conformal factors on the two vertices of the edge
	   *    \param e : edge e
	   */
	  virtual void _length( double u1, double u2, typename M::CEdge * e ) = 0;

	  /*! 
	   *	\brief Calculate corner angles using consine law
	   *    
	   */
	  virtual void _calculate_corner_angle();

	/*!
	 *	\brief hyperbolic cosine law
	 *  \param a,b,c are edge lengths, 
	 *  \return angle C
	 */
	double _cosine_law( double a, double b, double c ) ;

	/*!
	 *	\brief Cosine law, has to be defined in the derivated classes
	 *  \param a,b are edge lengths, c is the angle at C
	 *  \return edge length c
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
    virtual void _calculate_edge_weight() = 0;


	/*!
	 *	Set the target curvature on each vertex
	 */
    virtual void    _set_target_curvature();
	/*!
	 *	Curvature flow 
	 */
    virtual bool   _flow( double );

	 /*!
	  *	Newton's method to optimize the entropy energy
	  * \param threshold err bound
	  * \param step_length step length
	  */
     virtual void   _Newton( double threshold, double step_length );

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
CBaseHyperbolicRicciFlow<M>::CBaseHyperbolicRicciFlow( M * pMesh ):CBaseRicciFlow( pMesh)
{
};


//Calculate corner angle
template<typename M>
double CBaseHyperbolicRicciFlow<M>::_cosine_law( double a, double b, double c )
{
	double C;
	C = acos( (cosh(a) * cosh(b)-cosh(c))/(sinh(a)*sinh(b)) );
	return C;
};

//given edge lengths a,b and angle C, return edge length c

template<typename M>
double CBaseHyperbolicRicciFlow<M>::_inverse_cosine_law( double a, double b, double C )
{
  double c;
  c = cosh(a) * cosh(b) - sinh(a)*sinh(b) * cos(C) ;
  c = log( c+sqrt(c*c-1) );
  return c;
}

// u = log tanh (r/2)
template<typename M>
double CBaseHyperbolicRicciFlow<M>::_u2r( double u )
{  
    double e = exp(u);
    double r;
    r = log(  (1+e)/(1 - e ) );
    return r;
}

// u = log tanh (r/2)
template<typename M>
double CBaseHyperbolicRicciFlow<M>::_r2u( double r )
{  
    return log(tanh(r/2));
}



//set target curvature

template<typename M>
void CBaseHyperbolicRicciFlow<M>::_set_target_curvature()
{
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
    v->target_k() = 0;
  }
};


//gradient flow method

template<typename  M>
bool CBaseHyperbolicRicciFlow<M>::_flow( double error_threshold )
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

template<typename  M>
void CBaseHyperbolicRicciFlow<M>::_normalization( Eigen::VectorXd & x, int num )
{
	//for hyperbolic Ricci flow, do nothing
};




template<typename M>
void CBaseHyperbolicRicciFlow<M>::_calculate_corner_angle()
{
	
	for ( M::MeshFaceIterator fiter( m_pMesh); ! fiter.end(); fiter ++ )
  {
	  M::CFace * f = *fiter;

	  M::CHalfEdge * he[3];

      he[0] = m_pMesh->faceMostCcwHalfEdge( f ); 

	  for( int i = 0; i < 3; i ++ )
      {
        he[(i+1)%3] = m_pMesh->	faceNextCcwHalfEdge(he[i]);
      }

      double l[3];
      for(int i = 0; i < 3; i ++ )
      {
		  M::CEdge * e = m_pMesh->halfedgeEdge( he[i] );
          l[i] = e->length();
      }

      for(int i = 0; i < 3; i ++ )
      {
          he[(i+1)%3]->angle() = _cosine_law( l[(i+1)%3] ,l[(i+2)%3] ,l[i] );
      }
  }
}


//Newton's method for optimizing entropy energy

template<typename M>
void CBaseHyperbolicRicciFlow<M>::_Newton( double threshold, double step_length )
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
			v->huv()[0] = v->u();	
			v->u() = v->target_k() - v->k();
		}

		CPoisson<M> P(m_pMesh );
		//P.hyperbolic_solve2();
		P.hyperbolic_solve();
	

		for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
		{
			M::CVertex * v = *viter;
			v->u() = v->huv()[0] + v->u() * step_length;
		}
  }

};


//compute metric

template<typename M>
void CBaseHyperbolicRicciFlow<M>::_calculate_metric( double error_threshold, double step_length )
{
	_set_target_curvature();
	_Newton( error_threshold, step_length );

};






} //namespace RicciFlow
} //namespace MeshLib

#endif  _BASE_HYPERBOLIC_RICCI_FLOW_H_