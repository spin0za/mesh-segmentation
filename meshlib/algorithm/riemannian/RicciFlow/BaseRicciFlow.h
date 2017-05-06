/*! \file BaseRicciFlow.h
 *  \brief Base class for Ricci flow algorithm
 *  \author David Gu
 *  \date   documented on 10/17/2010
 *
 *	Algorithm for general Ricci Flow
 */

#ifndef _BASE_RICCI_FLOW_H_
#define _BASE_RICCI_FLOW_H_

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
#include "Operator/Operator.h"
#include "Laplace/Poisson.h"



namespace MeshLib
{

namespace RicciFlow
{

/*! \brief BaseClass CBaseRicciFlow
*
*	Algorithm for computing general Ricci flow
*/
template<typename M>
  class CBaseRicciFlow
  {
  public:
    /*! \brief CBaseRicciFlow constructor
	 *  \param pMesh the input mesh
	 */
    CBaseRicciFlow( M * pMesh );
    /*! \brief CBaseRicciFlow destructor
	 */
    ~CBaseRicciFlow();
	/*!	Computing the metric
	 */
	virtual void _calculate_metric( double threshold = 5e-4, double step_length = 1.0 );

  protected:
    /*!
     *	the input mesh
	 */
    M	  * m_pMesh;
	/*!
	 *	boundary of the input mesh
	 */
	typename M::CBoundary		  m_boundary;

  protected:
	  /*!
	   *	Calculate each edge length, has to be defined in the derivated classes
	   */
	  virtual void _length( double u1, double u2, typename M::CEdge * e )=0;
	  /*!
	   *	Calculate each edge length
	   */
   virtual void _calculate_edge_length();

	/*!
	 *	Cosine law, has to be defined in the derivated classes
	 */
	virtual double _cosine_law( double a, double b, double c ) = 0;
	/*!
	 *	Calculate corner angle
	 */
    virtual void _calculate_corner_angle();
	/*!
	 *	Calculate vertex curvature
	 */
    virtual void _calculate_vertex_curvature();
	/*!
	 *	Calculate vertex curvature error
	 */
   virtual double   _calculate_curvature_error();

	/*!
	 *	Calculate the edge weight
	 */
    virtual void _calculate_edge_weight() = 0;

	/*!
	 *	Set the target curvature on each vertex
	 */
    virtual void  _set_target_curvature() =  0;
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
	  *	Normalization
	  * \param du the du vector
	  * \param n dimension of the du vector
	  */
	 //virtual void _normalization( Eigen::VectorXd & du, int n ) = 0;
	 /*!
	  *	calculate hessian matrix Hessain 
	  * \param SparseMatrix
	  */
	 //virtual void _calculate_Hessain( Eigen::SparseMatrix<double> & pMatrix );
  };

//Constructor
template<typename M>
CBaseRicciFlow<M>::CBaseRicciFlow( M * pMesh ): m_pMesh( pMesh), m_boundary( pMesh )
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->u() = 0;
	}
};

//Destructor
template<typename M>
CBaseRicciFlow<M>::~CBaseRicciFlow( )
{
};

//Compute the edge length
template<typename M>
void CBaseRicciFlow<M>::_calculate_edge_length( )
{

	for( M::MeshEdgeIterator eiter(m_pMesh);  !eiter.end();  eiter++ )
  {
	  M::CEdge * e = *eiter;

	  M::CVertex * v1 = m_pMesh->edgeVertex1( e );
	  M::CVertex * v2 = m_pMesh->edgeVertex2( e );

      double u1 = v1->u();
      double u2 = v2->u();
	  
	  _length( u1, u2, e );
  }

};


//Calculate corner angle
template<typename M>
void CBaseRicciFlow<M>::_calculate_corner_angle( )
{
	COperator<M> pS(m_pMesh );
	pS._metric_2_angle();
};

//Calculate vertex curvature
template<typename M>
void CBaseRicciFlow<M>::_calculate_vertex_curvature()
{
	COperator<M> pS(m_pMesh );
	pS._corner_angle_2_vertex_curvature();
};


//compute curvature error

template<typename M>
double CBaseRicciFlow<M>::_calculate_curvature_error()
{
  double max_error = -1;
  M::CVertex * vert = NULL;

  for(M::MeshVertexIterator viter( m_pMesh); !viter.end() ; viter ++ )
   {
	   M::CVertex * v = *viter;
	   double k = v->target_k() - v->k();      
       k = fabs( k );
	   if( k > max_error )
	   {
		max_error = k;
		vert = v;
	   }
   }
   printf("Vertex id is %d\n", vert->id() );
   return max_error; 
};

//compute metric

template<typename M>
void CBaseRicciFlow<M>::_calculate_metric( double threshold, double step_length )
{

   _set_target_curvature();
   _Newton( threshold, step_length );     

};

//gradient flow method

template<typename M>
bool CBaseRicciFlow<M>::_flow( double error_threshold )
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
			double dif = v->target_k() - v->k();
            v->u() += dif * 2e-2;
		  }
    }
    return false;
};

//Newton's method for optimizing entropy energy

template<typename M>
void CBaseRicciFlow<M>::_Newton( double threshold, double step_length )
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
		P.solve2();
	

		for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
		{
			M::CVertex * v = *viter;
			v->u() = v->huv()[0] + v->u() * step_length;
		}

  }

};


//Newton's method for optimizing entropy energy
/*
template<typename V, typename E, typename F, typename H>
void CBaseRicciFlow<V,E,F,H>::_calculate_Hessain( Eigen::SparseMatrix<double> & M )
{
	std::vector<Eigen::Triplet<double> > M_coefficients;

	//set A
	for( CRicciFlowMesh<V,E,F,H>::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++  )
	{
	  E * e = *eiter;
	  V * v1 = m_pMesh->edgeVertex1( e );
	  V * v2 = m_pMesh->edgeVertex2( e );

	  M_coefficients.push_back( Eigen::Triplet<double>( v1->idx(), v2->idx(), -e->weight()) );
	  M_coefficients.push_back( Eigen::Triplet<double>( v2->idx(), v1->idx(), -e->weight()) );

	}
  
	for( CRicciFlowMesh<V,E,F,H>::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
	{
	    V * v = *viter;
		int idx = v->idx();
		double w = 0;
		for( CRicciFlowMesh<V,E,F,H>::VertexEdgeIterator veiter( v ); !veiter.end(); veiter ++ )
		{
		  E * pE = *veiter;
		  w += pE->weight();
		}
		M_coefficients.push_back( Eigen::Triplet<double>( idx, idx, w ) );
  }

  M.setFromTriplets(M_coefficients.begin(), M_coefficients.end());
}
*/
} //namespace RicciFlow
} //namespace MeshLib

#endif  _BASE_RICCI_FLOW_H_