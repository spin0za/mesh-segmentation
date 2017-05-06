/*! \file QCYamabeFlow.h
 * \brief Quasi-Conformal Euclidean Yamabe flow
 *  \author David Gu
 *  \date   documented on 03/08/2011
 *
 *	Algorithm for Quasi-Conformal Yamabe Flow
 */

#ifndef _QC_YAMABE_FLOW_H_
#define _QC_YAMABE_FLOW_H_

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
#include "QCRicciFlowMesh.h"
#include "BaseRicciFlow.h"
#include "Geometry/Circle.h"

namespace MeshLib
{

namespace RicciFlow
{
	/*! \brief CQCYamabeFlow class
	 *  
	 *  Algorithm for Quasi-Conformal Euclidean Yamabe flow
	 */
	 
 template<typename M>
 class CQCYamabeFlow : public CBaseRicciFlow<M>
  {
  public:
	  /*! \brief CQCYamabeFlow constructor
	   *
	   *  call base class constructor 
	   */
	  CQCYamabeFlow( M * pMesh );
	  /*! override the virtual function calculate_metric */
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
    void    _set_target_curvature();

  };

template<typename M>
CQCYamabeFlow<M>::CQCYamabeFlow( M * pMesh ): CBaseRicciFlow<M>( pMesh)
{

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		M::CVertex * v1 = m_pMesh->edgeVertex1(e);
		M::CVertex * v2 = m_pMesh->edgeVertex2(e);

	  CPoint2 p1 = v1->huv();
	  CPoint2 p2 = v2->huv();

	  std::complex<double> z1(p1[0], p1[1]);
	  std::complex<double> z2(p2[0], p2[1]);

	  std::complex<double> dz = z2 - z1;

	  std::complex<double> mu1 = v1->mu();
	  std::complex<double> mu2 = v2->mu();
	  std::complex<double> mu  = (mu1 + mu2)/2.0;

	  e->mu_length() = std::abs( dz + mu * std::conj( dz ));
	}


}

 //calculate edge length
template<typename M>
void CQCYamabeFlow<M>::_length( double u1, double u2, typename M::CEdge * e )
{
	  double l  = e->mu_length();
	  e->length() = exp(u1 ) * exp( u2 ) * l;
};

 //calculate edge weight

template<typename M>
void CQCYamabeFlow<M>::_calculate_edge_weight()
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
void CQCYamabeFlow<M>::_calculate_metric()
{

	  _set_target_curvature();
	  double error_threshold = 1e-6;
	  //double step_length = 0.5;
	  double step_length = 0.5;
	  _Newton( error_threshold, step_length );
		  
};      


//Euclidean cosine law
template<typename M>
double CQCYamabeFlow<M>::_cosine_law( double a, double b, double c )
{
          double cs =  ( a*a + b * b  - c * c )/( 2.0 * a * b );
          assert( cs <= 1.0 && cs >= -1.0 );
          return acos( cs );    
};

//normalization

template<typename M>
void CQCYamabeFlow<M>::_normalization( Eigen::VectorXd & x, int num )
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
void CQCYamabeFlow<M>::_set_target_curvature()
{
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
	  v->target_k() = 0;
  }
};

} //namespace RicciFlow
} //namespace MeshLib

#endif  _QC_YAMABE_FLOW_H_