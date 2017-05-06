#ifndef _POISSON_H_
#define _POISSON_H_

#include <map>
#include <vector>
#include <complex>

#include "Mesh/BaseMesh.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "mesh/boundary.h"
#include "Parser/parser.h"
#include <Eigen/Sparse>


namespace MeshLib
{
	template<typename M>
	class CPoisson
	{
		public:

		/*! \brief COperator constructor
		 *  \param pMesh the input mesh
		 */
		  CPoisson( M * pMesh ) { m_pMesh = pMesh; };
		/*! \brief COperator destructor
		 */
		  ~CPoisson(){};
		
		  /*!	\brief Solve Laplace Equation \Delta  u = b
		   *    \param input : the edge weight is set, the vertex u is set, which equals to b, 
		   *    \param output: the vertex u is set
		   */
		  void solve();
		  void solve2();
		  //for graph embedding
		  void solve3();
		  void hyperbolic_solve();
		  //void hyperbolic_solve2();

		protected:
			M * m_pMesh;
	};

template<typename M>
void CPoisson<M>::solve()
{
    
    int vid = 0;

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
    {
		M::CVertex * v = *viter;
        v->idx() = vid ++;
    }

	int num = m_pMesh->numVertices();

	std::vector<Eigen::Triplet<double> > A_coefficients;
	
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );

		int id1 = v1->idx();
		int id2 = v2->idx();

		double w = e->weight();
		
		A_coefficients.push_back( Eigen::Triplet<double>(id1,id2,-w) );
		A_coefficients.push_back( Eigen::Triplet<double>(id2,id1,-w) );

  }

	Eigen::VectorXd b( num );

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
    {
	   M::CVertex * v = *viter;
       int  id = v->idx();
	   b(id) = v->u();
	}

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
      double sum_w = 0;
	  M::CVertex * pV = *viter;
	  int id = pV->idx();

      for( M::VertexEdgeIterator veiter( pV ); !veiter.end();  ++veiter  )
      {
		  M::CEdge   * e = *veiter;
		  sum_w += e->weight();
       }

	  A_coefficients.push_back( Eigen::Triplet<double>(id,id,sum_w));
	}
    
	Eigen::SparseMatrix<double> A( num, num );
	A.setZero();
	A.setFromTriplets(A_coefficients.begin(), A_coefficients.end());


	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	std::cerr << "Eigen Decomposition" << std::endl;
	solver.compute(A);
	std::cerr << "Eigen Decomposition Finished" << std::endl;
	
	if( solver.info() != Eigen::Success )
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}

	Eigen::VectorXd x = solver.solve(b);
	if( solver.info() != Eigen::Success )
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		int  id = v->idx();
		v->u() = x(id);
	}
}


template<typename M>
void CPoisson<M>::solve2()
{
    
    int vid = 0;

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
    {
		M::CVertex * v = *viter;
        v->idx() = vid ++;
    }

	int num = m_pMesh->numVertices();

	std::vector<Eigen::Triplet<double> > A_coefficients;
	
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );

		int id1 = v1->idx();
		int id2 = v2->idx();

		double w = e->weight();
		
		A_coefficients.push_back( Eigen::Triplet<double>(id1,id2,-w) );
		A_coefficients.push_back( Eigen::Triplet<double>(id2,id1,-w) );

  }

	Eigen::VectorXd b( num );

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
    {
	   M::CVertex * v = *viter;
       int  id = v->idx();
	   b(id) = v->u();
	}

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
      double sum_w = 0;
	  M::CVertex * pV = *viter;
	  int id = pV->idx();

      for( M::VertexEdgeIterator veiter( pV ); !veiter.end();  ++veiter  )
      {
		  M::CEdge   * e = *veiter;
		  sum_w += e->weight();
       }

	  A_coefficients.push_back( Eigen::Triplet<double>(id,id,sum_w));
	}
    
	Eigen::SparseMatrix<double> A( num, num );
	A.setZero();
	A.setFromTriplets(A_coefficients.begin(), A_coefficients.end());


	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	std::cerr << "Eigen Decomposition" << std::endl;
	solver.compute(A);
	std::cerr << "Eigen Decomposition Finished" << std::endl;
	
	if( solver.info() != Eigen::Success )
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}

	Eigen::VectorXd x = solver.solve(b);
	if( solver.info() != Eigen::Success )
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		int  id = v->idx();
		v->u() = x(id);
	}

	//normalize

	double sum = 0;
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
	{
		M::CVertex * v = *viter;
		sum += v->u();
	}
	sum /= num;
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
	{
		M::CVertex * v = *viter;
		v->u()-= sum;
	}

};


/*
//for tangential hyperbolic ricci flow
template<typename M>
void CPoisson<M>::hyperbolic_solve()
{
 
    int vid = 0;
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
    {
		M::CVertex * v = *viter;
        v->idx() = vid ++;
    }

	int num = m_pMesh->numVertices();
	Eigen::VectorXd b( num );

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
    {
	   M::CVertex * v = *viter;
       int  id = v->idx();
	   b(id) = v->u();
	}


	std::vector<Eigen::Triplet<double> > _coefficients;

  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
  {
	  M::CVertex * v = *viter;

      double a_ii = 0;

	  for( M::VertexOutHalfedgeIterator hiter( m_pMesh, v ); !hiter.end(); ++ hiter )
	  {
		  M::CHalfEdge * he = *hiter;
	  	a_ii += he->c_w_ii();
	  }
	  _coefficients.push_back( Eigen::Triplet<double>( v->idx(), v->idx(), a_ii ) );
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

	  _coefficients.push_back( Eigen::Triplet<double>( v0->idx(), v1->idx(), a_ij ) );
	  _coefficients.push_back( Eigen::Triplet<double>( v1->idx(), v0->idx(), a_ij ) );
  }




	Eigen::SparseMatrix<double> A( num, num );
	A.setZero();
	A.setFromTriplets(_coefficients.begin(), _coefficients.end());


	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	std::cerr << "Eigen Decomposition" << std::endl;
	solver.compute(A);
	std::cerr << "Eigen Decomposition Finished" << std::endl;
	
	if( solver.info() != Eigen::Success )
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}

	Eigen::VectorXd x = solver.solve(b);
	if( solver.info() != Eigen::Success )
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		int  id = v->idx();
		v->u() = x(id);
	}
};
*/


//for tangential hyperbolic ricci flow
template<typename M>
void CPoisson<M>::hyperbolic_solve()
{
 
    int vid = 0;
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
    {
		M::CVertex * v = *viter;
        v->idx() = vid ++;
    }

	int num = m_pMesh->numVertices();
	Eigen::VectorXd b( num );

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
    {
	   M::CVertex * v = *viter;
       int  id = v->idx();
	   b(id) = v->u();
	}

	std::vector<Eigen::Triplet<double> > _coefficients;

  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
  {
	  M::CVertex * v = *viter;

      double a_ii = 0;

	  for( M::VertexInHalfedgeIterator hiter( m_pMesh, v ); !hiter.end(); ++ hiter )
	  {
		  M::CHalfEdge * he = *hiter;
	  	  a_ii += he->c_w_ii();
	  }
	  _coefficients.push_back( Eigen::Triplet<double> ( v->idx(), v->idx(), a_ii ) );
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

	  _coefficients.push_back( Eigen::Triplet<double> ( v0->idx(), v1->idx(), a_ij ) );
	  _coefficients.push_back( Eigen::Triplet<double> ( v1->idx(), v0->idx(), a_ij ) );

  }

	Eigen::SparseMatrix<double> A( num, num );
	A.setZero();
	A.setFromTriplets(_coefficients.begin(), _coefficients.end());


	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	std::cerr << "Eigen Decomposition" << std::endl;
	solver.compute(A);
	std::cerr << "Eigen Decomposition Finished" << std::endl;
	
	if( solver.info() != Eigen::Success )
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}

	Eigen::VectorXd x = solver.solve(b);
	if( solver.info() != Eigen::Success )
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		int  id = v->idx();
		v->u() = x(id);
	}
};



template<typename M>
void CPoisson<M>::solve3()
{
    
    int vid = 0;
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
    {
		M::CVertex * v = *viter;
		if( v->type() == 1 ) continue;
        v->idx() = vid ++;
    }

	int num = vid;

	std::vector<Eigen::Triplet<double> > A_coefficients;
	
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );
		
		if( v1->type() == 1 || v2->type() == 1 ) continue;

		int id1 = v1->idx();
		int id2 = v2->idx();

		double w = e->weight();
		
		A_coefficients.push_back( Eigen::Triplet<double>(id1,id2,-w) );
		A_coefficients.push_back( Eigen::Triplet<double>(id2,id1,-w) );

  }

	Eigen::VectorXd b( num );

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
    {
	   M::CVertex * v = *viter;
	   if( v->type() == 1 ) continue;
       int  id = v->idx();
	   b(id) = v->u();
	}

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
      double sum_w = 0;
	  M::CVertex * pV = *viter;
	  if( pV->type() == 1 ) continue;

	  int id = pV->idx();

      for( M::VertexEdgeIterator veiter( pV ); !veiter.end();  ++veiter  )
      {
		  M::CEdge   * e = *veiter;
		  M::CVertex * pV1 = m_pMesh->edgeVertex1( e );
		  M::CVertex * pV2 = m_pMesh->edgeVertex2( e );
		
		  if( pV1->type() == 1 ) continue;
		  if( pV2->type() == 1 ) continue;

		  sum_w += e->weight();
       }

	  A_coefficients.push_back( Eigen::Triplet<double>(id,id,sum_w));
	}
    
	Eigen::SparseMatrix<double> A( num, num );
	A.setZero();
	A.setFromTriplets(A_coefficients.begin(), A_coefficients.end());


	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	std::cerr << "Eigen Decomposition" << std::endl;
	solver.compute(A);
	std::cerr << "Eigen Decomposition Finished" << std::endl;
	
	if( solver.info() != Eigen::Success )
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}

	Eigen::VectorXd x = solver.solve(b);
	if( solver.info() != Eigen::Success )
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		if( v->type() == 1 ) continue;
		int  id = v->idx();
		v->u() = x(id);
	}

	//normalize

	double sum = 0;
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
	{
		M::CVertex * v = *viter;
		if( v->type() == 1 ) continue;
		sum += v->u();
	}
	sum /= num;
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
	{
		M::CVertex * v = *viter;
		if( v->type()==1 ) continue;
		v->u()-= sum;
	}

};


};
#endif _POISSON_H_