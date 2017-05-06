/*! \file DelaunayRefineMesh.h
*   \brief Mesh Refinement using Delaunay Triangulation
*   \date  Documented on 10/15/2010
*
*   Adding Steiner points to an existing mesh to improve the mesh quality. The input mesh must not have
*   flipped faces. In practice, we prefer to use harmonic mapping result as the input mesh, instead of Riemann mapping result.
*/

#ifndef  _DELAUNAY_REFINE_MESH_H_
#define  _DELAUNAY_REFINE_MESH_H_

#include <map>
#include <vector>
#include <queue>

#include "DelaunayMesh.h"

namespace MeshLib
{
  
/*! \brief CDeluanayRefineMesh class
*   
*   Mesh Refinement using Delaunay Triangulation
*/
template<typename V, typename E, typename F, typename H>
class CDelaunayRefineMesh : public CDelaunayMesh<V,E,F,H>
{
public:
	
	/*!
	* Refine a input mesh with uv coordinates
	* \param input_mesh input mesh file name
	*/
	void Refine( CDelaunayRefineMesh<V,E,F,H> & input_mesh );


};


//refine the mesh triangulation

template<typename V, typename E, typename F, typename H>
void CDelaunayRefineMesh<V,E,F,H>::Refine( CDelaunayRefineMesh<V,E,F,H> & input_mesh )
{
	_initialize( true ); //grading false = consider area; grading true, don't care area
	m_max_number_vertices = input_mesh.numVertices() * 2;

	CBoundary bnd( &input_mesh );
	std::vector<CLoop*> & loops = bnd.loops();

	CLoop * pL = loops[0];
	std::list<H*> & hs = pL->halfedges();

	std::vector<V*> segment_verts;
	
	for( std::list<H*>::iterator hiter = hs.begin() ; hiter != hs.end(); hiter ++ )
	{

		H * pH = *hiter;
		V * pV = halfedgeSource( pH );
		V * pW = NULL;
		_insert_vertex( pV->uv(), pW );
		pW->father() = pV->id();
		pW->point() = pV->point();
		segment_verts.push_back( pW );
		m_inputVertexNum ++;
	}
	for( size_t i = 0; i <  segment_verts.size(); i ++ )
	{
		V * v0 = segment_verts[i];
		V * v1 = segment_verts[(i+1)%segment_verts.size()];

		CSegment<V> * pS = _create_segment( v0, v1);
	}

	int c = 0;

	//input interior vertices
	for( MeshVertexIterator viter( &input_mesh ); !viter.end(); viter ++ )
	{
		V * pV = *viter;
		V * pW = NULL;
		
		if( pV->boundary()  ) continue;
		printf( "%d/%d\n", ++c, input_mesh.numVertices() );
		_insert_vertex( pV->uv(), pW );
		pW->father() = pV->id();
		pW->point() = pV->point();
	}

	//insert points
	_insert_missing_segments();
	_freeze_edges_on_segments();	
	_classify_inside_outside();

	RemoveBadTriangles( 30 );

}



typedef CDelaunayRefineMesh<CDelaunayVertex, CDelaunayEdge, CDelaunayFace, CDelaunayHalfEdge> CDRTMesh;


}
#endif  _DELAUNAY_REFINE_MESH_