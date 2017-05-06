/*!
*      \file PunctureMesh.h
*      \brief Mesh for punching a hole in the center
*	   \author David Gu
*      \date Documented on 10/12/2010
*
*/
/*******************************************************************************
*      Puncture Mesh
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       Puncture a hole in the center of the mesh
* 
*       David Gu June 27, 2008,  gu@cs.stonybrook.edu
*
*******************************************************************************/

#ifndef  _PUNCTURE_MESH_H_
#define  _PUNCTURE_MESH_H_

#include "Mesh/Vertex.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "Mesh/HalfEdge.h"

#include "Mesh/BaseMesh.h"
#include "Mesh/iterators.h"
#include "Mesh/boundary.h"

namespace MeshLib
{
namespace Topology
{
/*! \brief CPunctureVertex class
*
*	Vertex class for punching a hole in the center
*   trait: whether the vertex has been touched m_touched
*/
  class CPunctureVertex : public  CVertex
  {
  protected:
    bool        m_touched;
  public:
	  bool    &   touched() { return m_touched; };
  };

/*---------------------------------------------------------------------------------------------------------------------------------------

	Puncture Mesh

----------------------------------------------------------------------------------------------------------------------------------------*/
/*!	\brief CPunctureMesh class
*
*	Mesh class for punching a hole in the center
*/
template<typename V, typename E, typename F, typename H>
class CPunctureMesh : public CBaseMesh<V,E,F,H>
{
public:
	typedef V CVertex;
	typedef E CEdge;
	typedef F CFace;
	typedef H CHalfEdge;

	typedef CBoundary<V,E,F,H> CBoundary;
	typedef CLoop<V,E,F,H> CLoop;
	typedef MeshFaceIterator<V,E,F,H> MeshFaceIterator;
	typedef MeshVertexIterator<V,E,F,H> MeshVertexIterator;
	typedef MeshEdgeIterator<V,E,F,H> MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
	typedef VertexFaceIterator<V,E,F,H> VertexFaceIterator;
public:
};


typedef CPunctureMesh<CPunctureVertex, CEdge, CFace, CHalfEdge> CPMesh;

} //namespace Topology
} //namespace MeshLib
#endif  _PUNCTURE_MESH_H_