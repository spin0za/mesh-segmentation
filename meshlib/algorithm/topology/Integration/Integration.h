/*! \file Integration.h
*   \brief Algorithm for integrating holomorphic 1-forms on a mesh
*	\author David Gu
*   \date Documented on 10/12/2010
*/

/********************************************************************************************************************************
*
*      Integration  Class
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       Integration Class
* 
*       David Gu June 27, 2008,  gu@cs.stonybrook.edu
*
*
*      Input:
*         
*           Holomorphic 1-form mesh, and a foundamental domain mesh ( a topological disk )
*
*      Output:
*
*           Integrate the holomorphic 1-form on the fundamental domain, to compute the vertex uv coordinates
*
*********************************************************************************************************************************/


/*---------------------------------------------------------------------------------------------------------------------------------
#include <math.h>
#include "Integration/Integration.h"


using namespace MeshLib;

int main( int argc, char * argv[] )
{
	CIMesh holo_mesh;
	CIMesh fund_mesh;

	holo_mesh.read_m( argv[1] );
	fund_mesh.read_m( argv[2] );

	CIntegration integrator( & holo_mesh, & fund_mesh );
	integrator._integrate();

	fund_mesh.write_m( argv[3] );

	return 0;
}
----------------------------------------------------------------------------------------------------------------------------------*/
#ifndef _INTEGRATION_H_
#define _INTEGRATION_H_

#include <queue>
#include "IntegrationMesh.h"

namespace MeshLib
{

namespace Topology
{
/*! \brief CIntegration class
*
*	Integrates a holomorphic 1-form on a mesh, and get he uv coordiantes for each vertex
*/ 
template<typename M>
class CIntegration
{
public:
	/*! CIntegration constructor
	* \param pForm input holomorphic 1-form 
	* \param pDomain output the integration results are stored in the vertex uv field
	*/
	CIntegration( M * pForm, M * pDomain );
	/*! CIntegration destructor */
	~CIntegration();
	/*! Integrate holomorphic 1-form on the mesh */
	void _integrate();

private:
	/*! Holomorphic 1-form */
	M * m_pForm;
	/*! The integration domain mesh, the integration results are stored in the vertex uv field */
	M * m_pDomain;

};

//CIntegration constructor
//\param pForm the input holomorphic 1-form
//\param pDomain the domain mesh with vertex uv
template<typename M>
CIntegration<M>::CIntegration( M * pForm, M * pDomain )
{
	m_pForm = pForm;
	m_pDomain = pDomain;

};

//CIntegration destructor
template<typename M>
CIntegration<M>::~CIntegration()
{
};


//Compute the integration
//Using breadth first search 
template<typename M>
void CIntegration<M>::_integrate()
{
	M::CVertex * head = NULL;

	for( M::MeshVertexIterator viter( m_pDomain ); !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->touched() = false;
		head = v;
	}


	std::queue<M::CVertex*> vqueue;
	head->uv() = CPoint2(0,0);
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
				e->duv() = we->duv();
			}
			else
			{
				e->duv() = CPoint2(0,0) - we->duv();
			}

			if( tail == v2 )
			{
				tail->uv()  = head->uv() + e->duv();
			}
			else
			{
				tail->uv()  = head->uv() - e->duv();
			}
		}
	}
	
}



} //namespace Topology
} //namespace MeshLib
#endif