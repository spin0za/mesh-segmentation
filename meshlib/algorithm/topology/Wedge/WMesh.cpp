/*!
*      \file WMesh.cpp
*      \brief Implement CWMesh class
*	   \author David Gu
*      \date Document 10/11/2010
*
*		Slice a mesh along its sharp edges to form a new open mesh.
*
*/

#include "WMesh.h"

using namespace MeshLib::Topology;


/*! CWedge constructor
 *	\param half_edge the first corner of the wedge
 */

//////////////////////////////////////////////////////////////////////////////////////////
//	constructor
//////////////////////////////////////////////////////////////////////////////////////////
CWedge::CWedge(CSMesh * pMesh, CWedgeHalfEdge * half_edge)
	  :m_mesh(pMesh),m_vertex( m_mesh->halfedgeVertex( half_edge ) )
{
	CWedgeEdge * edge = m_mesh->halfedgeEdge( half_edge );

	//if the half_edge isBoundary, then it must be the most Ccwone

	if( edge->boundary() || edge->sharp()  )
	{
		if( edge->boundary() )
		{
			CWedgeVertex * pW = m_mesh->halfedgeTarget( half_edge );
			assert( m_mesh->edgeHalfedge( edge,0 ) == m_mesh->vertexMostCcwInHalfEdge( pW ) );
		}

		/*
		 *	build a half wedge
		 */

		CWedgeHalfEdge * he = m_mesh->edgeHalfedge( edge, 0 );
		
		if( he->vertex() != m_vertex )
			he = m_mesh->halfedgeSym( he );

		m_mostCcwHEdge = he;
		CWedgeHalfEdge * he_prev;

		do{
			he_prev = he;
			he = m_mesh->vertexNextClwInHalfEdge( he );
		}while( he != NULL && ! m_mesh->halfedgeEdge( he )->sharp() );
		
		m_mostClwHEdge = he_prev;
	}
	else
	{
		/*
		 *	build a full wedge
		 */

		CWedgeHalfEdge * he = m_mesh->vertexMostCcwInHalfEdge( m_vertex );

		m_mostCcwHEdge = he;
		m_mostClwHEdge = m_mesh->vertexNextCcwInHalfEdge( he );
	}

};


/*-----------------------------------------------------------------------------------------------------------------

	Constructor/ Destructor

------------------------------------------------------------------------------------------------------------------*/
/*! CWMesh constructor 
*/

CWMesh::CWMesh()
{
}

/*! CWMesh constructor 
*  \param pMesh the input mesh
*/

CWMesh::CWMesh( CSMesh * pMesh )
{
	m_pMesh = pMesh;
};

/*! CWMesh destructor 
*  release all the memories for the list of wedges.
*/


CWMesh::~CWMesh()
{
  for( std::list<CWedge*>::iterator  witer =  m_wedges.begin() ; witer != m_wedges.end(); witer ++ )
	{
		CWedge * w = *witer;
		delete w;
	}

};

/*-----------------------------------------------------------------------------------------------------------------

	Slice the mesh along sharp edges to form a new open mesh

------------------------------------------------------------------------------------------------------------------*/
/*! Slice the mesh along sharp edges to form a new open mesh.
 */

void CWMesh::Slice()
{
	_construct();
	_convert();
};


/*-----------------------------------------------------------------------------------------------------------------

	Compute the topological valence of a vertex

------------------------------------------------------------------------------------------------------------------*/

int CWMesh::__topovalence( CWedgeVertex * vertex )
{

	int valence = 0;

	for( CSMesh::VertexEdgeIterator veiter( vertex ); !veiter.end(); ++ veiter)
	{
		CWedgeEdge * edge = *veiter;
		if( edge->boundary() || edge->sharp() )
		{
			valence ++;
		}
	}

	return valence;
};

/*-----------------------------------------------------------------------------------------------------------------

	Construct a mesh, each vertex is an wedge

------------------------------------------------------------------------------------------------------------------*/

void CWMesh::_construct()
{
	
	for( CSMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		CWedgeVertex * vertex = *viter;
		vertex->valence() = __topovalence( vertex );

		if(  vertex->valence() == 0 )
		{
			//find most ccw HEdge of vertex
			CWedgeHalfEdge * he = m_pMesh->vertexMostCcwInHalfEdge(vertex);
			assert( he != NULL );

			//generate a new wedge
			CWedge * wedge = new CWedge( m_pMesh, he );
			assert( wedge );
			m_wedges.push_back( wedge );

			for( CSMesh::VertexInHalfedgeIterator vhiter( m_pMesh, vertex ); !vhiter.end(); vhiter ++ )
			{
				CWedgeHalfEdge * pH = *vhiter;
				pH->wedge() = wedge;
			}
		}
		else
		{

			//boundary vertex

			//assume this rotate Clw, if the vertex is on the boundary, 
			//the first edge is the most Ccw Edge

			CWedgeHalfEdge * he = m_pMesh->vertexMostCcwInHalfEdge( vertex );
			CWedge		   * head = NULL;
			CWedgeEdge     * pE = m_pMesh->halfedgeEdge( he );

			while( !pE->boundary() && !pE->sharp() )
			{
				  he = m_pMesh->vertexNextClwInHalfEdge( he );
				  pE = m_pMesh->halfedgeEdge( he );
			}

			CWedgeHalfEdge * anchor = he;

			while( true )
			{

				if( m_pMesh->halfedgeEdge( he )->boundary() || m_pMesh->halfedgeEdge( he )->sharp() )
				{
					//generate a wedge
					CWedge * wedge = new CWedge( m_pMesh, he );
					assert( wedge );
					m_wedges.push_back( wedge );

					//assign all corners' wedge
					for( he = wedge->mostCcwHalfEdge(); he != wedge->mostClwHalfEdge(); he = m_pMesh->vertexNextClwInHalfEdge( he ) )
					{
						he->wedge() = wedge;
						pE = m_pMesh->halfedgeEdge( he );
					}

					assert( he == wedge->mostClwHalfEdge() );
					he->wedge() = wedge;
				}

				he = m_pMesh->vertexNextClwInHalfEdge( he ); //he->clw_rotate_about_target();
				if( he == NULL && m_pMesh->isBoundary( vertex ) ) break;
				if( he == anchor ) break;
			}

		}


	}

};

/*-----------------------------------------------------------------------------------------------------------------

	Convert wedge mesh to common mesh

------------------------------------------------------------------------------------------------------------------*/

void CWMesh::_convert()
{

	int ind = 1;

	for( std::list<CWedge*>::iterator witer =  m_wedges.begin( ); witer != m_wedges.end();  witer++ )
	{
		CWedge * wedge = * witer;
		CWedgeVertex * wvertex = m_wmesh.createVertex( ind ++ );
		assert( wvertex );

		wvertex->string()= wedge->vertex()->string();
		wvertex->point() = wedge->vertex()->point();
		wedge->wvertex() = wvertex;
		wvertex->father() = wedge->vertex()->id();
		//For Computing Tight Fundamental Group Generators
		//wvertex->parent() = (wedge->vertex()->parent()!= 0)?wedge->vertex()->parent():wedge->vertex()->id();
	}

	int find = 1;
	for( CSMesh::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); ++ fiter )
	{
		CFace * f = *fiter;
		std::vector<CWedgeVertex*>    v;

		for( CSMesh::FaceHalfedgeIterator fhiter( f ); !fhiter.end(); ++ fhiter )
		{
			CWedgeHalfEdge * he = *fhiter;
			CWedge *    w  = he->wedge();
			assert( w );
			v.push_back( w->wvertex() );
		}

		CFace *wf = m_wmesh.createFace(v,find++);
		wf->string() = f->string();
	}

	for( CSMesh::MeshEdgeIterator eiter( &m_wmesh); !eiter.end(); ++ eiter )
	{
		CWedgeEdge * e = *eiter;
		
		CWedgeVertex * v1 = m_wmesh.edgeVertex1( e );
		CWedgeVertex * v2 = m_wmesh.edgeVertex2( e );
		
		int f1 = v1->father();
		int f2 = v2->father();
		
		CWedgeVertex * fv1 = m_pMesh->idVertex( f1 );
		CWedgeVertex * fv2 = m_pMesh->idVertex( f2 );
		
		CWedgeEdge   * pFE = m_pMesh->vertexEdge( fv1, fv2 );
		e->string() = pFE->string();

		if( e->boundary() )
		{
			v1->boundary() = true;
			v2->boundary() = true;
		}
	}


	//Arrange the boundary half_edge of boundary vertices, to make its halfedge
	//to be the most ccw in half_edge

	for(std::list<CWedgeVertex*>::iterator viter = m_wmesh.vertices().begin();  viter != m_wmesh.vertices().end() ; ++ viter )
	{
		CWedgeVertex *     v = *viter;
		if( !v->boundary() ) continue;

		CWedgeHalfEdge * he = m_wmesh.vertexHalfedge( v );

		while( m_wmesh.halfedgeSym( he ) != NULL )
		{
			he = m_wmesh.vertexNextCcwInHalfEdge( he );
		}
		v->halfedge() = he;
	}

	//copy corner information
	for( CSMesh::MeshFaceIterator fiter( &m_wmesh ); !fiter.end(); ++ fiter )
	{
		CFace * f = *fiter;
		for( CSMesh::FaceHalfedgeIterator fhiter( f ); !fhiter.end(); ++ fhiter )
		{
			CWedgeHalfEdge * he = *fhiter;
			
			CWedgeVertex *   target = m_wmesh.halfedgeTarget( he );
			CWedgeVertex *   source = m_wmesh.halfedgeSource( he );
			
			int  ft = target->father();
			int  fs = source->father();

			CWedgeVertex *   pWT = m_pMesh->idVertex( ft );
			CWedgeVertex *   pWS = m_pMesh->idVertex( fs );

			CWedgeHalfEdge * pWH = m_pMesh->vertexHalfedge( pWS, pWT );
			if( pWH == NULL )
			{
				printf("Face ID %d Error\n", f->id() );
			}
			assert( pWH != NULL );
			he->string() = pWH->string();
		}
	}


};

