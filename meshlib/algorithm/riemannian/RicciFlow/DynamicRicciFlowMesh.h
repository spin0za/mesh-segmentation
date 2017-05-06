/*!
*      \file  DynamicRicciFlowMesh.h
*      \brief Dynamic Mesh for Ricci Flow
*	   \author David Gu
*      \date Documented 02/08/2013
*
*/

#ifndef  _DYNAMIC_RICCI_FLOW_MESH_H_
#define  _DYNAMIC_RICCI_FLOW_MESH_H_

#include <map>
#include <vector>

#include "Mesh/BaseMesh.h"
#include "Mesh/DynamicMesh.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "mesh/boundary.h"
#include "Parser/parser.h"
#include "Parser/traits_io.h"
#include "RicciFlowMesh.h"
#include "Mesh/DynamicDelaunayMesh.h"

namespace MeshLib
{
namespace RicciFlow
{

/*-------------------------------------------------------------------------------------------------------------------------------------

	Ricci flow Dynamic Mesh Class

--------------------------------------------------------------------------------------------------------------------------------------*/
/*!
 *	\brief CDynamicRicciFlowMesh class
 *
 *	Mesh class for dynamic Ricci Flow mesh purpose
 */
template<typename V, typename E, typename F, typename H>
class CDynamicRicciFlowMesh : public CDynamicDelaunayMesh<V,E,F,H>
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
	typedef VertexFaceIterator<V,E,F,H> VertexFaceIterator;
	typedef VertexInHalfedgeIterator<V,E,F,H> VertexInHalfedgeIterator;
	typedef VertexOutHalfedgeIterator<V,E,F,H> VertexOutHalfedgeIterator;
	typedef FaceHalfedgeIterator<V,E,F,H> FaceHalfedgeIterator;
	typedef FaceEdgeIterator<V,E,F,H> FaceEdgeIterator;
	typedef MeshFaceIterator<V,E,F,H> MeshFaceIterator;
	typedef FaceVertexIterator<V,E,F,H> FaceVertexIterator;
	typedef MeshHalfEdgeIterator<V,E,F,H> MeshHalfEdgeIterator;
public:
};


typedef CDynamicRicciFlowMesh<CRicciFlowVertex, CRicciFlowEdge, CRicciFlowFace, CRicciFlowHalfEdge> CDRFMesh;	

} //namespace RicciFlow
} //namespace MeshLib
#endif  _DYNAMIC_RICCI_FLOW_MESH_H_