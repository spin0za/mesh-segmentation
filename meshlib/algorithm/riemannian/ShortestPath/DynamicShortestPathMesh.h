/*!
*      \file DynamicShortestPathMesh.h
*      \brief Dynamic Mesh for Shortest Path
*	   \author David Gu
*      \date Document 10/12/2013
*
*/

#ifndef  _DYNAMIC_SHORTEST_PATH_MESH_H_
#define  _DYNAMIC_SHORTEST_PATH_MESH_H_

#include <map>
#include <vector>

#include "Mesh/BaseMesh.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "mesh/boundary.h"
#include "Parser/traits_io.h"
#include "ShortestPathMesh.h"

namespace MeshLib
{
/*!	\brief CDynamicShortestPathMesh class
 *
 *	Dynamic Mesh class for computing shortest paths between two boundary loops
 */

template<typename V, typename E, typename F, typename H>
class CDynamicShortestPathMesh : public CDynamicMesh<V,E,F,H>
{
public:
	typedef V CVertex;
	typedef E CEdge;
	typedef F CFace;
	typedef H CHalfEdge;

	typedef CBoundary<V,E,F,H> CBoundary;
	typedef CLoop<V,E,F,H> CLoop;
	typedef MeshVertexIterator<V,E,F,H> MeshVertexIterator;
	typedef MeshEdgeIterator<V,E,F,H> MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
public:
};


typedef CDynamicShortestPathMesh<CSPVertex, CSPEdge, CFace, CHalfEdge> CDSPMesh;

};
#endif  _DYNAMIC_SHORTEST_PATH_MESH_H_