/*!
*      \file DynamicWMesh.h
*      \brief Algorithm for slicing a mesh along sharp edges to form a new mesh
*	   \author David Gu
*      \date Document 10/11/2010
*
*		Slice a mesh along sharp edges to form a new mesh. Input is a mesh with sharp edges,
*		the output is the open mesh sliced along the sharp edges. 
*/





/*------------------------------------------------------------------------------------------------------------------------------------

//Example Code

#include "mesh/mesh.h"
#include "Wedge/WMesh.h"

using namespace MeshLib;

int main( int argc, char * argv[] )
{
	CMesh mesh;
	mesh.read_m( argv[1] );

	CWMesh wmesh( & mesh );
    wmesh.Slice();

    wmesh.wmesh()->write_m( argv[2] );

  return 0;
}

------------------------------------------------------------------------------------------------------------------------------------*/

#ifndef  _DYNAMIC_WMesh_H_
#define  _DYNAMIC_WMesh_H_


#include <assert.h>
#include <math.h>

#include <iostream>
#include <map>
#include <list>

#include "DynamicWedgeMesh.h"


namespace MeshLib{

namespace Topology
{

class CBaseWedge
{
};

/*!
 *  \brief CDynamicWedge class
 *
 *   A wedge is a union of corners sharing the same apex vertex. The corners attaching to the same vertex are
 *   partitioned to different wedges by sharp edges or boundary edges.
 */
template <typename M>
class CDynamicWedge  : public CBaseWedge
{
public:
	/*! CDynamicWedge constructure */
	CDynamicWedge(){};
	/*! CDynamicWedge destructor */
	~CDynamicWedge(){};

	/*! CDynamicWedge constructor
	*	\param pMesh the input mesh
	*	\param half_edge the first corner in this wedge
	*/
	CDynamicWedge( M * pMesh, typename M::CHalfEdge *half_edge);//constructor 
	/*!	The vertex of the current wedge
	*/
	typename M::CVertex *  vertex() { return m_vertex; };						//the vertex of the wedge
	/*!	The most clockwise corner of the current wedge
	 */
	typename M::CHalfEdge * mostClwHalfEdge(){ return m_mostClwHEdge;};			//return the most Clw HEdge of this wedge
	/*!	The most counter clockwise corner of the current wedge
	 */
	typename M::CHalfEdge * mostCcwHalfEdge(){return m_mostCcwHEdge;};			//return the most Ccw HEdge of this wedge
	/*!	Each wedge corresponds to a vertex in wmesh, the associated wmesh vertex of current wedge
	 */
	typename M::CVertex * & wvertex() { return m_wvertex; };					//the new vertex on the wedge solid

protected:
	/*!	The current wedge mesh
	*/
	M					* m_mesh;
	/*!	vertex of all the corners in the wedge
	*/
	typename M::CVertex			* m_vertex;
	/*!	The most CCW corner of the wedge.
	 */
	typename M::CHalfEdge		* m_mostCcwHEdge;
	/*! The most CLW corner of the wedge.
	 */
	typename M::CHalfEdge        * m_mostClwHEdge;
	/*!	New vertex on the wedge mesh, corresponding to this wedge.
	 */
	typename M::CVertex			* m_wvertex;	//new vertex on wedge solid

};

//Wedge Solid, Support mesh slicing
/*!
 *	\brief CWMesh class
 *
 *  Converting a mesh with sharp edges to a new mesh, such that 
 *  1. each wedge becomes a new vertex
 *  2. each face of the old mesh becomes a new face in the following way: each corner of the old face belongs to an wedge,
 *     the three wedges are connected to a new face.
 */
template<typename M>
class CDynamicWMesh : public M
{

public:
	/*! CDynamicWMesh constructor */
	//constructor and destructor
	CDynamicWMesh(){};
	/*! CDynamicWMesh constructor
	* \param pMesh the input mesh with sharp edges
	*/
	CDynamicWMesh( M * pMesh ){ m_pMesh = pMesh; };
	/*! CWMesh destructor */
	~CDynamicWMesh();
	/*! The newly constructed mesh */
	M * wmesh() { return &m_wmesh; };	
	/*! Slice the input mesh along the sharp edges. */
	void   Slice();

private:
	/*! Construct a halfedge structure, each vertex is a wedge. */
	void		_construct();
	/*! Convert the halfedge structure with wedge vertices to a common mesh. */
	void		_convert();
	/*! The topological valence of the vertex, number of sharp edges or boundary edges,
	 *  \param vertex input vertex
	 */
	int		    __topovalence( typename M::CVertex * vertex );			//compute topological valence of vertex

	/*! The input mesh. */
	M *			m_pMesh;
	/*! The output converted mesh. */
	M		    m_wmesh;
	/*! list of wedges. */
	std::list<CDynamicWedge<M>*> m_wedges;							//buffer for wedges
};

template <typename M>
//////////////////////////////////////////////////////////////////////////////////////////
//	constructor
//////////////////////////////////////////////////////////////////////////////////////////
CDynamicWedge<M>::CDynamicWedge(M * pMesh, typename M::CHalfEdge * half_edge)
	  :m_mesh(pMesh),m_vertex( m_mesh->halfedgeVertex( half_edge ) )
{
	M::CEdge * edge = m_mesh->halfedgeEdge( half_edge );

	//if the half_edge isBoundary, then it must be the most Ccwone

	if( edge->boundary() || edge->sharp()  )
	{
		if( edge->boundary() )
		{
			M::CVertex * pW = m_mesh->halfedgeTarget( half_edge );
			assert( m_mesh->edgeHalfedge( edge,0 ) == m_mesh->vertexMostCcwInHalfEdge( pW ) );
		}

		/*
		 *	build a half wedge
		 */

		M::CHalfEdge * he = m_mesh->edgeHalfedge( edge, 0 );
		
		if( he->vertex() != m_vertex )
			he = m_mesh->halfedgeSym( he );

		m_mostCcwHEdge = he;
		M::CHalfEdge * he_prev;

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

		M::CHalfEdge * he = m_mesh->vertexMostCcwInHalfEdge( m_vertex );

		m_mostCcwHEdge = he;
		m_mostClwHEdge = m_mesh->vertexNextCcwInHalfEdge( he );
	}

};


/*! CDynamicWMesh destructor 
*  release all the memories for the list of wedges.
*/

template<typename M>
CDynamicWMesh<M>::~CDynamicWMesh()
{
  for( std::list<CDynamicWedge<M>*>::iterator  witer =  m_wedges.begin() ; witer != m_wedges.end(); witer ++ )
	{
		CDynamicWedge<M> * w = *witer;
		delete w;
	}

};


/*-----------------------------------------------------------------------------------------------------------------

	Slice the mesh along sharp edges to form a new open mesh

------------------------------------------------------------------------------------------------------------------*/
/*! Slice the mesh along sharp edges to form a new open mesh.
 */
template<typename M>
void CDynamicWMesh<M>::Slice()
{
	_construct();
	_convert();
};

/*-----------------------------------------------------------------------------------------------------------------

	Compute the topological valence of a vertex

------------------------------------------------------------------------------------------------------------------*/
template<typename M>
int CDynamicWMesh<M>::__topovalence( typename M::CVertex * vertex )
{

	int valence = 0;

	for( M::VertexEdgeIterator veiter( vertex ); !veiter.end(); ++ veiter)
	{
		M::CEdge * edge = *veiter;
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
template<typename M>
void CDynamicWMesh<M>::_construct()
{
	
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * vertex = *viter;
		vertex->valence() = __topovalence( vertex );

		if(  vertex->valence() == 0 )
		{
			//find most ccw HEdge of vertex
			M::CHalfEdge * he = m_pMesh->vertexMostCcwInHalfEdge(vertex);
			assert( he != NULL );

			//generate a new wedge
			CDynamicWedge<M> * wedge = new CDynamicWedge<M>( m_pMesh, he );
			assert( wedge );
			m_wedges.push_back( wedge );

			for( M::VertexInHalfedgeIterator vhiter( m_pMesh, vertex ); !vhiter.end(); vhiter ++ )
			{
				M::CHalfEdge * pH = *vhiter;
				pH->wedge() = wedge;
			}
		}
		else
		{

			//boundary vertex

			//assume this rotate Clw, if the vertex is on the boundary, 
			//the first edge is the most Ccw Edge

			M::CHalfEdge * he = m_pMesh->vertexMostCcwInHalfEdge( vertex );
			CDynamicWedge<M>		   * head = NULL;
			M::CEdge     * pE = m_pMesh->halfedgeEdge( he );

			while( !pE->boundary() && !pE->sharp() )
			{
				  he = m_pMesh->vertexNextClwInHalfEdge( he );
				  pE = m_pMesh->halfedgeEdge( he );
			}

			M::CHalfEdge * anchor = he;

			while( true )
			{

				if( m_pMesh->halfedgeEdge( he )->boundary() || m_pMesh->halfedgeEdge( he )->sharp() )
				{
					//generate a wedge
					CDynamicWedge<M> * wedge = new CDynamicWedge<M>( m_pMesh, he );
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
template<typename M>
void CDynamicWMesh<M>::_convert()
{

	int ind = 1;
	//generate new vertices
	for( std::list<CDynamicWedge<M>*>::iterator witer =  m_wedges.begin( ); witer != m_wedges.end();  witer++ )
	{
		CDynamicWedge<M> * wedge = * witer;
		M::CVertex * wvertex = new M::CVertex();
		assert( wvertex );
		m_verts.push_back( wvertex );
		wvertex->id() = ind ++;
		m_map_vert.insert( std::pair<int,M::CVertex*>( wvertex->id(), wvertex) );
		wvertex->halfedge() = NULL;
		wvertex->string()= wedge->vertex()->string();
		wvertex->point() = wedge->vertex()->point();
		wvertex->father() = wedge->vertex()->id();
		wedge->wvertex() = wvertex;
	}

	ind = 1;

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); ++ eiter )
	{
		M::CEdge * e = *eiter;
		M::CHalfEdge * ph[2];
		M::CHalfEdge * ch[2];

		//generate half edges attaching to the current edge
		
		for( int k = 0; k < 2; k ++ )
		{
			ph[k] = edgeHalfedge(e, k);
			if( ph[k] != NULL )
			{
				ch[k] = new M::CHalfEdge();
				assert( ch[k] );
				ph[k]->child() = ch[k];
				
				//connecting halfedge -> vertex pointer
				CDynamicWedge<M>* wedge = (CDynamicWedge<M>*) ph[k]->wedge();
				ch[k]->vertex() = wedge->wvertex();
				wedge->wvertex()->halfedge() = ch[k];
				ch[k]->string() = ph[k]->string();
			}
			else
			{
				ch[k] = NULL;
			}
		}

		//connecting halfedge edge
		if( !e->sharp() )
		{
			M::CEdge * ce = new M::CEdge();
			assert( ce );
			ce->id() = ind ++;

			m_edges.push_back( ce );
			for( int k = 0; k < 2; k ++ )
			{
				ce->halfedge(k) = ch[k];
				if( ch[k] != NULL )
					ch[k]->edge() = ce;
			}
			
			//std::cout << "Edge id " << ce->id() << std::endl;

			ce->string() = e->string();
		}

		if( e->sharp() )
		{

			M::CEdge * ce[2];

			for( int k = 0; k < 2; k ++ )
			{
				ce[k]  = new M::CEdge();
				assert( ce[k] );
				ce[k]->id() = ind ++;
				m_edges.push_back( ce[k] );
				ce[k]->halfedge(0) = ch[k];
				ce[k]->halfedge(1) = NULL;
				ch[k]->edge() = ce[k];
				ce[k]->string() = e->string();
				
				//std::cout << "Sharp Edge id " << ce[k]->id() << std::endl;

			}
		}
	}

/*
	for( std::list<M::CEdge*>::iterator eiter = m_edges.begin(); eiter != m_edges.end(); ++ eiter )
	{
		M::CEdge * e = *eiter;
		std::cout << "Edge " << e->id() << std::endl;
		
		M::CHalfEdge * h = edgeHalfedge( e, 0 );
		M::CVertex   * v = halfedgeVertex( h );
		std::cout << "Target " << v->id() << std::endl;

		h = edgeHalfedge( e, 1 );
		if( h == NULL ) continue;

		v = halfedgeVertex( h );
		std::cout << "Target " << v->id() << std::endl;
	}
*/
	int find = 1;
	for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); ++ fiter )
	{
		M::CFace * f = *fiter;
		std::vector<M::CHalfEdge*>    phs;

		for( M::FaceHalfedgeIterator fhiter( f ); !fhiter.end(); ++ fhiter )
		{
			M::CHalfEdge * he = *fhiter;
			phs.push_back( he->child() );
		}

		M::CFace *wf = new M::CFace();
		assert( wf );
		wf->id() = find ++;
		wf->string() = f->string();
		m_faces.push_back( wf );
		m_map_face.insert( std::pair<int,M::CFace*>( wf->id(), wf ) );

		for( size_t k = 0; k < 3; k ++ )
		{
			phs[k]->face() = wf;
			wf->halfedge() = phs[k];
			phs[k]->he_next() = phs[(k+1)%3];
			phs[k]->he_prev() = phs[(k+2)%3];
		}
	}
/* Debug purpose only
	for( std::list<M::CEdge*>::iterator eiter = m_edges.begin(); eiter != m_edges.end(); ++ eiter )
	{
		M::CEdge * e = *eiter;
		std::cout << "Edge " << e->id() << std::endl;
		
		M::CHalfEdge * h = edgeHalfedge( e, 0 );
		M::CVertex   * v = halfedgeVertex( h );
		std::cout << "Target " << v->id() << std::endl;

		h = edgeHalfedge( e, 1 );
		if( h == NULL ) continue;

		v = halfedgeVertex( h );
		std::cout << "Target " << v->id() << std::endl;

		M::CVertex * v1 = edgeVertex1( e );
		M::CVertex * v2 = edgeVertex2( e );

		std::cout << "Vertex 1 " << v1->id() << std::endl;
		std::cout << "Vertex 2 " << v2->id() << std::endl;
	}
*/
	for( std::list<M::CVertex*>::iterator viter = m_verts.begin(); viter != m_verts.end(); ++ viter )
	{
		M::CVertex * v = *viter;
		v->boundary() = false;
	}

	for( std::list<M::CEdge*>::iterator eiter = m_edges.begin(); eiter != m_edges.end(); ++ eiter )
	{
		M::CEdge * e = *eiter;
		//std::cout << e->id() << std::endl;

		if( e->boundary() )
		{
			M::CVertex * v1 = edgeVertex1( e );
			M::CVertex * v2 = edgeVertex2( e );
			v1->boundary() = true;
			v2->boundary() = true;
		}
	}

	//Arrange the boundary half_edge of boundary vertices, to make its halfedge
	//to be the most ccw in half_edge

	for(std::list<M::CVertex*>::iterator viter = m_verts.begin();  viter != m_verts.end() ; ++ viter )
	{
		M::CVertex *     v = *viter;
		if( !v->boundary() ) continue;

		M::CHalfEdge * he = vertexHalfedge( v );

		while( halfedgeSym( he ) != NULL )
		{
			he = vertexNextCcwInHalfEdge( he );
		}
		v->halfedge() = he;
	}

};

}//name space Topology
}//name space MeshLib

#endif //_MESHLIB_SOLID_H_ defined



