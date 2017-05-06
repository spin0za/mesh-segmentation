/*!
*      \file BaseHolomorphicForm.h
*      \brief Algorithm for computing holomorphic 1-forms, Hodge star operator
*	   \author David Gu
*      \date Document 10/12/2010
*
*		Algorithm for computing holomorphic 1-forms, Hodge Star operator
*
*/
/********************************************************************************************************************************
*
*      BaseHolomorphic 1-form Class
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       Compute the holomorphic 1-forms
* 
*       David Gu June 27, 2008,  gu@cs.stonybrook.edu
*
*
*      Input:
*         
*           Original mesh, the mesh cut through a shortest path connecting one inner boundary and the exterior boundary
*
*      Output:
*
*           Closed non-exact Harmonic 1-form. and the mesh with UV coordinates.
*
*********************************************************************************************************************************/

/*---------------------------------------------------------------------------------------------------------------------------------
#include "HolomorphicForm/HolomorphicForm.h"

using namespace MeshLib;

int main( int argc, char * argv[] )
{
	if( strcmp( argv[1], "-holomorphic_form" ) == 0 )
	{
	std::list<CHoloFormMesh*> meshes;


	for( int i = 2; i < argc; i ++ )
	{
		CHoloFormMesh * pMesh = new CHoloFormMesh;
		assert( pMesh );
		pMesh->read_m( argv[i] );
		meshes.push_back( pMesh );
	}

	CHolomorphicForm form( meshes );
	form.conjugate();

	int id = 0;
	for( std::list<CHoloFormMesh*>::iterator miter = meshes.begin(); miter != meshes.end(); miter++)
	{
		CHoloFormMesh * pM = *miter;
		std::stringstream iss;
		iss << "holomorphic_form_" << id++ << ".m";
		pM->write_m( iss.str().c_str() );
	}

	for( std::list<CHoloFormMesh*>::iterator miter = meshes.begin(); miter != meshes.end(); miter++)
	{
		CHoloFormMesh * pM = *miter;
		delete pM;
	}
	
}
----------------------------------------------------------------------------------------------------------------------------------*/
#ifndef _BASE_HOLOMORPHIC_FORM_H_
#define _BASE_HOLOMORPHIC_FORM_H_

#include  <math.h>
#include <queue>
#include <list>
#include <vector>
#include "HolomorphicFormMesh.h"
#include "Structure/Structure.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>


namespace MeshLib
{

namespace Holomorphy
{

/*! \brief CWedgeOperator class
 *
 *	Wedge Star Operator
 */
template<typename M>
class CWedgeOperator  
{
public:
	/*! CWedgeOperator constructor
	 * \param solid0 harmonic 1-form \f$\omega_0\f$
	 * \param solid1 harmonic 1-form \f$\omega_1\f$
	 */
	CWedgeOperator( M * mesh0, M * mesh1);
	/*! CWedgeOperator destructor
	*/
	~CWedgeOperator();
	/*! wedge product 
	 * \return 
	 *	\f[
				\int \omega_0 \wedge \omega_1
        \f]
	 */
	double wedge_product();
	/*! wedge product 
	 * \return 
	 *	\f[
				\int \omega_0 \wedge {}^*\omega_1
        \f]
	 */
	double wedge_star_product();
private:
	/*! Two input harmonic 1-forms $\f\omega_0,\omega_1\f$. */
	M * m_pMesh[2];
};

/*! \brief CHolomorphicForm class
 *
 *	Compute holomorphic forms on a mesh
 */
template<typename M>
class CBaseHolomorphicForm
{
public:
	/*! CHolomorphicForm constructor
	*	\param meshes the list of meshes with stores the harmonic 1-forms, which form the basis 
	*   of the first cohomology group of the mesh
	*/
	CBaseHolomorphicForm( std::list<M*> & meshes );
	/*! CHolomorphicForm destructor
	*/
	~CBaseHolomorphicForm();
	/*!	The list of meshes storing harmonic form bases
	*/
	std::vector<M*> & meshes() { return m_meshes; };

	/*! Compute the conjugate harmonic 1-forms for the base harmonic 1-forms by solving the equation,
		Assume \f$ {}^*\omega_i = \sum_j \lambda_{ij} \omega_j \f$, then
		\f[
			\int_M \omega_k \wedge {}^*\omega_i = \sum_j \lambda_{ij} \int_M \omega_k \wedge \omega_j
		\f]
		the left hand side can be computed using face-vector represenation.
	*/
	void conjugate();

protected:
	/*!	The list of meshes storing harmonic form bases
	*/
	std::vector<M*> m_meshes;
	/*! Convert harmonic 1-form from edge-function to face-vector representation 
	 */
	//void convert( CHoloFormMesh * pMesh, CHoloFormFace * f );

	virtual void _angle_structure() = 0;
};

template<typename M>
double CWedgeOperator<M>::wedge_product()
{

	double p = 0;
	for( M::MeshFaceIterator fiter( m_pMesh[0] ); !fiter.end(); ++ fiter )
	{
		M::CFace * f0 = *fiter;
		M::CFace * f1 = m_pMesh[1]->idFace( f0->id() );
		
		std::vector<M::CHalfEdge*> h0,h1;
		for( M::FaceHalfedgeIterator fhiter( f0 ); !fhiter.end(); ++fhiter )
		{
			M::CHalfEdge * h = *fhiter;
			h0.push_back( h );
		}

		for( M::FaceHalfedgeIterator fhiter( f1 ); !fhiter.end(); ++fhiter )
		{
			M::CHalfEdge * h = *fhiter;
			h1.push_back( h );
		}

		M::CVertex *  s0 = m_pMesh[0]->halfedgeSource( h0.front() );
		for( size_t i = 0; i < 3; i ++ )
		{
			M::CHalfEdge * h = h1.front();
			M::CVertex *   s = m_pMesh[1]->halfedgeSource( h );
			if( s->id() != s0->id() )
			{
				h1.erase( h1.begin() );
				h1.push_back( h );
			}
			else
			{
				break;
			}
		}

		std::vector<double> du0,du1;
		for( size_t i = 0; i < 3; i ++ )
		{
			M::CHalfEdge * h = h0[i];
			M::CEdge * e = m_pMesh[0]->halfedgeEdge( h );
			double du = ( h == m_pMesh[0]->edgeHalfedge( e, 0 )  )? e->du():-e->du();
			du0.push_back( du );	
		}

		for( size_t i = 0; i < 3; i ++ )
		{
			M::CHalfEdge * h = h1[i];
			M::CEdge * e = m_pMesh[1]->halfedgeEdge( h );
			double du = ( h == m_pMesh[1]->edgeHalfedge( e, 0 )  )? e->du():-e->du();
			du1.push_back( du );	
		}


		
		//	 |du0[0] du0[1] du0[2]|
		//	 |du1[0] du1[1] du1[2]| /6
		//	 |1      1      1     |
		 
		p += ( du0[1] * du1[2] - du0[2] * du1[1] ) /2.0;
	}
	return p;
};

//Hodge Star operator
template<typename M>
double CWedgeOperator<M>::wedge_star_product()
{

	double p = 0;
	for( M::MeshFaceIterator fiter( m_pMesh[0] ); !fiter.end(); ++ fiter )
	{
		M::CFace * f0 = *fiter;
		M::CFace * f1 = m_pMesh[1]->idFace( f0->id() );

		std::vector<M::CHalfEdge*> h0,h1;
		for( M::FaceHalfedgeIterator fhiter( f0 ); !fhiter.end(); ++fhiter )
		{
			M::CHalfEdge * h = *fhiter;
			h0.push_back( h );
		}

		for( M::FaceHalfedgeIterator fhiter( f1 ); !fhiter.end(); ++fhiter )
		{
			M::CHalfEdge * h = *fhiter;
			h1.push_back( h );
		}

		M::CVertex *  s0 = m_pMesh[0]->halfedgeSource( h0.front() );
		for( size_t i = 0; i < 3; i ++ )
		{
			M::CHalfEdge * h = h1.front();
			M::CVertex *   s = m_pMesh[1]->halfedgeSource( h );
			if( s->id() != s0->id() )
			{
				h1.erase( h1.begin() );
				h1.push_back( h );
			}
			else
			{
				break;
			}
		}

		std::vector<double> du0,du1,theta;
		for( size_t i = 0; i < 3; i ++ )
		{
			M::CHalfEdge * h = h0[i];
			M::CEdge * e = m_pMesh[0]->halfedgeEdge( h );
			double du = ( h == m_pMesh[0]->edgeHalfedge( e, 0 )  )? e->du():-e->du();
			du0.push_back( du );	
			theta.push_back( h0[(i+1)%3]->angle() );
		}

		for( size_t i = 0; i < 3; i ++ )
		{
			M::CHalfEdge * h = h1[i];
			M::CEdge * e = m_pMesh[1]->halfedgeEdge( h );
			double du = ( h == m_pMesh[1]->edgeHalfedge( e, 0 )  )? e->du():-e->du();
			du1.push_back( du );	
		}
	

		p += 0.5 * du0[0] * du1[0] * cos( theta[0] )/sin( theta[0] );
		p += 0.5 * du0[1] * du1[1] * cos( theta[1] )/sin( theta[1] );
		p += 0.5 * du0[2] * du1[2] * cos( theta[2] )/sin( theta[2] );
	}
	return p;
};

//CWedgeOperator constructor
//\param mesh0, mesh1 two input harmonic 1-forms
template<typename M>
CWedgeOperator<M>::CWedgeOperator( M * mesh0, M * mesh1)
{
	m_pMesh[0] = mesh0;
	m_pMesh[1] = mesh1;
};

//CWedgeOperator destructor
template<typename M>
CWedgeOperator<M>::~CWedgeOperator()
{
};

//CBaseHolomorphicForm constructor
//\param meshes are the basis of harmonic 1-form group
template<typename M>
CBaseHolomorphicForm<M>::CBaseHolomorphicForm( std::list<M*> & meshes )
{
	for( std::list<M*>::iterator miter = meshes.begin(); miter != meshes.end(); miter ++ )
	{
		M * pM = *miter;
		m_meshes.push_back( pM );
	}
};

//CBaseHolomorphicForm destructor
template<typename M>
CBaseHolomorphicForm<M>::~CBaseHolomorphicForm()
{
};


// Compute the conjugate harmonic 1-forms for the base harmonic 1-forms
template<typename M>
void CBaseHolomorphicForm<M>::conjugate()
{

	_angle_structure();

	int n = m_meshes.size();


	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);

	for(int i = 0;  i < n; i ++ )
	{
		for( int j = 0;  j < n; j ++ )
		{
				CWedgeOperator<M> wo( m_meshes[i], m_meshes[j] );
				A(i,j) = wo.wedge_product();
		}
	}



	for(int i = 0; i < n ; i ++ )
	{
		Eigen::VectorXd b(n);
		
		for( int j = 0; j < n ; j ++ )
		{
			CWedgeOperator<M> wo( m_meshes[i], m_meshes[j] );
			b(j) = wo.wedge_star_product();
		}

		
		//Eigen::VectorXd x = solver.solve(b);
		Eigen::VectorXd x = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

		
		std::cout << "Hodge Star Coefficients: ";
		for( int j = 0; j < n ; j ++ )
		{
			std::cout << x(j) << " ";
		}
		std::cout << std::endl;

		for( M::MeshEdgeIterator eiter( m_meshes[i] ); !eiter.end(); ++ eiter )
		{
			M::CEdge * e = *eiter;
			e->duv()[0] = e->du();
			e->duv()[1] = 0;
		
			int id1 = m_meshes[i]->edgeVertex1(e)->id();
			int id2 = m_meshes[i]->edgeVertex2(e)->id();

			for( int k = 0; k < n ; k ++ )
			{
				M::CVertex * w1 = m_meshes[k]->idVertex( id1 );
				M::CVertex * w2 = m_meshes[k]->idVertex( id2 );

				M::CEdge * edge = m_meshes[k]->vertexEdge( w1, w2 );
				e->duv()[1] += edge->du() * x(k);
			}

		}

	}

};





} //namespace Holomorphy
} //namespace MeshLib
#endif