/*! \file DynamicBaseRicciFlow.h
 *  \brief Base class for Dynamic Ricci flow algorithm
 *  \author David Gu
 *  \date   documented on 02/08/2013
 *
 *	Algorithm for general Dynamic Ricci Flow (Ricci flow with edge swap, such that 
 *	the mesh is always Delaunay )
 */

#ifndef _DYNAMIC_BASE_RICCI_FLOW_H_
#define _DYNAMIC_BASE_RICCI_FLOW_H_

#include <map>
#include <vector>
#include <Eigen/Sparse>

#include "Mesh/BaseMesh.h"
#include "Mesh/DynamicMesh.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "mesh/boundary.h"
#include "Parser/parser.h"
#include "RicciFlowMesh.h"
#include "DynamicRicciFlowMesh.h"
#include "Operator/Operator.h"

namespace MeshLib
{

namespace RicciFlow
{

/*! \brief BaseClass CDynamicBaseRicciFlow
*
*	Algorithm for computing general Dynamic Ricci flow
*/
template<typename M>
  class CDynamicBaseRicciFlow
  {
  public:
    /*! \brief CDynamicBaseRicciFlow constructor
	 *  \param pMesh the input mesh
	 */
    CDynamicBaseRicciFlow( M * pMesh );
    /*! \brief CDynamicBaseRicciFlow destructor
	 */
    ~CDynamicBaseRicciFlow();
	/*!	Computing the metric
	 */
	virtual void _calculate_metric();


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
    void _calculate_edge_length();

	/*!
	 *	Cosine law, has to be defined in the derivated classes
	 */
	virtual double _cosine_law( double a, double b, double c ) = 0;
	/*!
	 *	Calculate corner angle
	 */
    void _calculate_corner_angle();
	/*!
	 *	Calculate vertex curvature
	 */
    void _calculate_vertex_curvature();
	/*!
	 *	Calculate vertex curvature error
	 */
   double   _calculate_curvature_error();

	/*!
	 *	Calculate the edge weight
	 */
    virtual void _calculate_edge_weight() = 0;

	/*!
	 *	Set the target curvature on each vertex
	 */
    virtual void  _set_target_curvature() =  0;
	 /*!
	  *	Dynamic Newton's method to optimize the entropy energy
	  * if the step length is too big, shrink it to a half
	  * \param threshold err bound
	  * \param step_length step length
	  */
     virtual void   _Dynamic_Newton( double threshold, double step_length );
	 /*!
	  *	Dynamic Newton's method to optimize the entropy energy one step
	  * if the step length is too big, shrink it to a half
	  * \param step_length step length
	  */
     virtual bool   _forward( double step_length );
	 virtual bool   _validate();

	 /*!
	  *	Normalization
	  * \param du the du vector
	  * \param n dimension of the du vector
	  */
	 virtual void _normalization( Eigen::VectorXd & du, int n ) = 0;
	 /*!
	  *	calculate hessian matrix Hessain 
	  * \param SparseMatrix
	  */
	 virtual void _calculate_Hessain( Eigen::SparseMatrix<double> & pMatrix );
	 /*!
	  *	preserve Delaunay Triangulation during the flow
	  *
	  */
	 virtual void _preserve_Delaunay() = 0;
  };

/*! 
	Constructor, 
	a. set vertex index
	b. set edge length, induced Euclidean metric
	c. set vertex conformal factor 0
*/

template<typename M>
CDynamicBaseRicciFlow<M>::CDynamicBaseRicciFlow( M * pMesh ): m_pMesh( pMesh), m_boundary( pMesh )
{
	int idx = 0;
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->idx() = idx ++; 
	}
	
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->u() = 0; 
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		e->length() = m_pMesh->edgeLength( e );
	}
	
};

//Destructor
template<typename M>
CDynamicBaseRicciFlow<M>::~CDynamicBaseRicciFlow( )
{
};

//Compute the edge length
template<typename M>
void CDynamicBaseRicciFlow<M>::_calculate_edge_length( )
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
void CDynamicBaseRicciFlow<M>::_calculate_corner_angle( )
{
	COperator<M> pS( m_pMesh );
	pS._metric_2_angle();
/*
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
		  if( l[(i+1)%3] + l[(i+2)%3] <= l[i] )
		  {
			  std::cerr <<"CDynamicBaseRicciFlow<M>::_calculate_corner_angle::breaks triangle inequality" << std::endl; 
		  }
          he[(i+1)%3]->angle() = _cosine_law( l[(i+1)%3] ,l[(i+2)%3] ,l[i] );
      }
  }
*/
};

//Calculate vertex curvature
template<typename M>
void CDynamicBaseRicciFlow<M>::_calculate_vertex_curvature()
{
	COperator<M> pS( m_pMesh );
	pS._angle_2_curvature();

/*
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
      double k  = (v->boundary() )? PI: PI * 2;
	  for( M::VertexInHalfedgeIterator vh( m_pMesh, v ); !vh.end();  ++vh )
      {
		  M::CHalfEdge * he = *vh;
          k -= he->angle();
      }
      v->k() = k;
  }
*/
};


//compute curvature error

template<typename M>
double CDynamicBaseRicciFlow<M>::_calculate_curvature_error()
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
void CDynamicBaseRicciFlow<M>::_calculate_metric()
{

   _set_target_curvature();
  //double error = 1e-6;
  double threshold = 5e-4;
  double step_length = 1.0;
   //_Newton( threshold, step_length );     

   _Dynamic_Newton( threshold, step_length );
};


//Newton's method for optimizing entropy energy

template<typename M>
void CDynamicBaseRicciFlow<M>::_calculate_Hessain( Eigen::SparseMatrix<double> & M )
{
	std::vector<Eigen::Triplet<double> > M_coefficients;

	//set A
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++  )
	{
		M::CEdge * e = *eiter;
		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );
	  //important
	  if( v1 == v2 ) continue;

	  if( v1->idx() == 0 || v2->idx() == 0 ) continue;

	  M_coefficients.push_back( Eigen::Triplet<double>( v1->idx()-1, v2->idx()-1, -e->weight()) );
	  M_coefficients.push_back( Eigen::Triplet<double>( v2->idx()-1, v1->idx()-1, -e->weight()) );

	}
  
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
	{
		M::CVertex * v = *viter;
		int idx = v->idx();
		if( idx == 0 ) continue;

		double w = 0;
		for( M::VertexEdgeIterator veiter( v ); !veiter.end(); veiter ++ )
		{
			M::CEdge * pE = *veiter;
			M::CVertex * v1 = m_pMesh->edgeVertex1( pE );
			M::CVertex * v2 = m_pMesh->edgeVertex2( pE );
		  if( v1 == v2 ) continue;

		  w += pE->weight();
		}
		M_coefficients.push_back( Eigen::Triplet<double>( idx-1, idx-1, w ) );
  }

  M.setFromTriplets(M_coefficients.begin(), M_coefficients.end());
}


//verify if the mesh is geometric, meaning if all the faces satisfy triangle inequality

template<typename M>
bool CDynamicBaseRicciFlow<M>::_validate()
{
 		_calculate_edge_length();
		if( !m_pMesh->_is_metric() )
		{
			std::cerr << "Invalid Metric" << std::endl;
			return false;
		}
		return true;
}

template<typename M>
bool CDynamicBaseRicciFlow<M>::_forward( double step_length )
{
 		_calculate_edge_length();
		_calculate_corner_angle();
		_calculate_vertex_curvature();
		_calculate_edge_weight();
		
		int num = m_pMesh->numVertices();
		Eigen::SparseMatrix<double>  M( num-1, num-1 );
		M.setZero();
		_calculate_Hessain( M );

		Eigen::VectorXd b(num-1);

		for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
		{
			M::CVertex * v = *viter;
			int idx = v->idx();
			if( idx == 0 ) continue;
			b(idx-1) = v->target_k() - v->k();

			for( M::VertexEdgeIterator veiter( v ); !veiter.end(); veiter ++ )
			{
				M::CEdge * pE = *veiter;
				M::CVertex * v1 = m_pMesh->edgeVertex1( pE );
				M::CVertex * v2 = m_pMesh->edgeVertex2( pE );
			  if( v1->idx() == 0 || v2->idx()==0 )
				b(idx-1) += pE->weight();
			}
		}



		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		std::cerr << "Eigen Decomposition" << std::endl;
		solver.compute(M);
		std::cerr << "Eigen Decomposition Finished" << std::endl;
	
		if( solver.info() != Eigen::Success )
		{
			std::cerr << "Waring: Eigen decomposition failed" << std::endl;
			return false;
		}

		Eigen::VectorXd y = solver.solve(b);
		if( solver.info() != Eigen::Success )
		{
			std::cerr << "Waring: Eigen decomposition failed" << std::endl;
			return false;
		}

		Eigen::VectorXd x(num);
		x(0) = 1.0;
		for( int i = 1; i < num; i ++ )
		{
			x(i) = y(i-1);
		}
		_normalization( x, num );

		Eigen::VectorXd old_u(num);
		for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
		{
			M::CVertex * v = *viter;
			int idx = v->idx();
			old_u(idx) = v->u();
		}

		do{
			for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
			{
				M::CVertex * v = *viter;
				int idx = v->idx();
				v->u() += x(idx) * step_length;
			}
			
			if( _validate() ) break;

			for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
			{
				M::CVertex * v = *viter;
				int idx = v->idx();
				v->u() = old_u(idx);
			}

			step_length /= 2.0;

		}while( true );

		return true;
}

//Newton's method for optimizing entropy energy

template<typename M>
void CDynamicBaseRicciFlow<M>::_Dynamic_Newton( double threshold, double step_length )
{
	
	_preserve_Delaunay();

  	while( true )
	{
		_forward( step_length );
		_preserve_Delaunay();

		double error =  _calculate_curvature_error();
		printf("Newton's Method: Current error is %f\r\n", error );
		if( error < threshold) break;
  }

};

} //namespace RicciFlow
} //namespace MeshLib

#endif  _DYNAMIC_BASE_RICCI_FLOW_H_