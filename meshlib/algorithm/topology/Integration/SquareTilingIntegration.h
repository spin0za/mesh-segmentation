/*! \file SquareTilingIntegration.h
*   \brief Algorithm for integrating a combinatorial harmonic 1-form on a mesh to obtain a square tiling
*	\author David Gu
*   \date Documented on 10/12/2010
*/

/********************************************************************************************************************************
*
*      Square Tiling Integration  Class
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       Square Tiling Integration Class
* 
*       David Gu October 25, 2014,  gu@cs.stonybrook.edu
*
*
*      Input:
*         
*           Combinatorial Harmonic 1-form mesh, and a foundamental domain mesh ( a topological disk )
*
*      Output:
*
*           Integrate a combinatorial harmonic 1-form on the fundamental domain, each edge is mapped to a square
*
*********************************************************************************************************************************/


#ifndef _SQUARE_TITLING_INTEGRATION_H_
#define _SQUARE_TITLING_INTEGRATION_H_

#include <queue>
#include "IntegrationMesh.h"

namespace MeshLib
{

namespace Topology
{
/*! \brief CSquareTilingIntegration class
*
*	Integrates a combinatorial 1-form on a mesh, and each edge is mapped to a rectangle
*/ 
template<typename M>
class CSquareTilingIntegration
{
public:
	/*! CSquareTilingIntegration constructor
	* \param pForm input holomorphic 1-form 
	* \param pDomain output the integration results are stored in the vertex uv field
	*/
	CSquareTilingIntegration( M * pForm, M * pDomain );
	/*! CIntegration destructor */
	~CSquareTilingIntegration();
	/*! Integrate holomorphic 1-form on the mesh */
	void _integrate();

private:
	/*! Holomorphic 1-form */
	M * m_pForm;
	/*! The integration domain mesh, the integration results are stored in the vertex uv field */
	M * m_pDomain;

	/*! integrate on the vertices */
	void _integrate_vertex();
	/*! integrate on the faces, the dua mesh vertices */
	void _integrate_face();
};

//CSquareTilingIntegration constructor
//\param pForm the input holomorphic 1-form
//\param pDomain the domain mesh with vertex uv
template<typename M>
CSquareTilingIntegration<M>::CSquareTilingIntegration( M * pForm, M * pDomain )
{
	m_pForm = pForm;
	m_pDomain = pDomain;

};

//CIntegration destructor
template<typename M>
CSquareTilingIntegration<M>::~CSquareTilingIntegration()
{
};


//Compute the integration
//Using breadth first search 
template<typename M>
void CSquareTilingIntegration<M>::_integrate_vertex()
{
	M::CVertex * head = NULL;

	for( M::MeshVertexIterator viter( m_pDomain ); !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->touched() = false;
		head = v;
	}


	std::queue<M::CVertex*> vqueue;
	head->u() = 0;
	head->touched()  = true;

	vqueue.push( head );

	while( !vqueue.empty() )
	{
		head = vqueue.front();
		vqueue.pop();

		for( M::VertexEdgeIterator veiter( head ); !veiter.end(); veiter ++ )
		{
			M::CEdge * e = *veiter;
			M::CVertex * v1 = m_pDomain->edgeVertex1( e );
			M::CVertex * v2 = m_pDomain->edgeVertex2( e );


			M::CVertex * tail = ( head != v1 )?v1:v2;
			if( tail->touched() ) continue;

			tail->touched() = true;
			vqueue.push( tail );

			int id1 = v1->father();
			//if there is no "father" field for the vertex, then directly use the vertex id 
			if( id1 == 0 ) id1 = v1->id();
			int id2 = v2->father();
			//if there is no "father" field for the vertex, then directly use the vertex id 
			if( id2 == 0 ) id2 = v2->id();

			M::CVertex * w1 = m_pForm->idVertex( id1 );
			M::CVertex * w2 = m_pForm->idVertex( id2 );

			M::CEdge * we = m_pForm->vertexEdge( w1, w2 );

			if( m_pForm->edgeVertex1( we ) == w1 )
			{
				e->du() = we->du();
			}
			else
			{
				e->du() = 0 - we->du();
			}

			if( tail == v2 )
			{
				tail->u()  = head->u() + e->du();
			}
			else
			{
				tail->u()  = head->u() - e->du();
			}
		}
	}

};

//Compute the integration
//Using breadth first search 
template<typename M>
void CSquareTilingIntegration<M>::_integrate_face()
{
	M::CFace * head = NULL;

	for( M::MeshFaceIterator fiter( m_pDomain ); !fiter.end(); fiter ++ )
	{
		M::CFace * f = *fiter;
		f->touched() = false;
		head = f;
	}


	std::queue<M::CFace*> fqueue;
	head->u() = 0;
	head->touched()  = true;

	fqueue.push( head );

	while( !fqueue.empty() )
	{
		head = fqueue.front();
		fqueue.pop();

		for( M::FaceHalfedgeIterator fhiter( head ); !fhiter.end(); fhiter ++ )
		{
			M::CHalfEdge * h = *fhiter;
			M::CEdge     * e =  m_pDomain->halfedgeEdge( h );
			M::CHalfEdge * sh = m_pDomain->halfedgeSym( h );
			if( sh == NULL ) continue;

			M::CFace     *  f = m_pDomain->halfedgeFace( sh );
			if( f->touched() ) continue;

			f->touched() = true;
			fqueue.push( f );

			M::CVertex * v1 = m_pDomain->halfedgeSource( h );
			M::CVertex * v2 = m_pDomain->halfedgeTarget( h );

			int id1 = v1->father();
			//if there is no "father" field for the vertex, then directly use the vertex id 
			if( id1 == 0 ) id1 = v1->id();
			int id2 = v2->father();
			//if there is no "father" field for the vertex, then directly use the vertex id 
			if( id2 == 0 ) id2 = v2->id();

			M::CVertex * w1 = m_pForm->idVertex( id1 );
			M::CVertex * w2 = m_pForm->idVertex( id2 );

			M::CEdge * we = m_pForm->vertexEdge( w1, w2 );

			double du = 0;

			if( m_pForm->edgeVertex1( we ) == w1 )
			{
				du = we->du();
			}
			else
			{
				du =  -we->du();
			}

			f->u() = head->u() - du;
		}
	}

};

//Compute the integration
//Using breadth first search 
template<typename M>
void CSquareTilingIntegration<M>::_integrate()
{
	_integrate_vertex();
	_integrate_face();
};


} //namespace Topology
} //namespace MeshLib
#endif