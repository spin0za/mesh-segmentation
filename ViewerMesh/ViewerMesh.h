/*! \file ViewerMesh.h
*   \brief Viewer Mesh
*   \author David Gu
*   \date   documented on 10/12/2010
*
*   Mesh for viewer 
*/
#ifndef  _VIEWER_MESH_H_
#define  _VIEWER_MESH_H_

#include <map>
#include <vector>

#include "Mesh/Vertex.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "Mesh/HalfEdge.h"

#include "Mesh/BaseMesh.h"
#include "Mesh/boundary.h"
#include "Mesh/iterators.h"
#include "Parser/parser.h"
#include "Parser/traits_io.h"

namespace MeshLib
{
/*! \brief CViewerFace class
*
*	Face class for viewer
*   Trait: face normal
*/
  class CViewerFace: public CFace
  {
  protected:
    /*! face normal */
    CPoint m_normal;
  public:
    /*! face normal */
	  CPoint & normal() { return m_normal; };
  };

  /*! \brief CViewerVertex class
  *   
  *   Vertex class for viewer
  *   Trait : vertex rgb color
  */
  class CViewerVertex : public  CVertex
  {
  protected:
	  /*! vertex rgb color */
    CPoint   m_rgb;
	
  public:
	  /*! CViewerVertex Constructor */
	  CViewerVertex()
	  {
		  m_rgb = CPoint(1,1,1); //default color is white
	  }
	  /*! vertex rgb color */
	  CPoint  & rgb()    { return m_rgb;    };
	
	  /*! read vertex rgb, uv from vertex string */
	  void _from_string();
  };
 
 // read vertex rgb, uv from vertex string 
 inline void CViewerVertex::_from_string()
 {
		  CParser parser( m_string );
		
		  for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
		  {
			  CToken * token = *iter;
			  if( token->m_key == "uv" )
			  {
				  token->m_value >> m_uv;
			  }
			  else
			  if( token->m_key == "rgb" )
			  {
				  token->m_value >> m_rgb;
			  }
		  }
  };



/*-------------------------------------------------------------------------------------------------------------------------------------

	Viewer Mesh

--------------------------------------------------------------------------------------------------------------------------------------*/
/*! \brief CViewerMesh class
*
*	mesh class for viewer
*
*/
 template<typename V, typename E, typename F, typename H>
class CViewerMesh : public CBaseMesh<V,E,F,H>
{
public:
	typedef CBoundary<V,E,F,H> CBoundary;
	typedef CLoop<V,E,F,H> CLoop;
	typedef MeshVertexIterator<V,E,F,H> MeshVertexIterator;
	typedef MeshFaceIterator<V,E,F,H>   MeshFaceIterator;
	typedef MeshEdgeIterator<V,E,F,H>   MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef FaceVertexIterator<V,E,F,H> FaceVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
	typedef VertexFaceIterator<V,E,F,H> VertexFaceIterator;
	
public:
};

typedef CViewerMesh<CViewerVertex, CEdge, CViewerFace, CHalfEdge> CVMesh;

unsigned long long CVMesh::m_input_traits  = VERTEX_UV| VERTEX_RGB;
unsigned long long CVMesh::m_output_traits = 0;

}
#endif  VIEWER_MESH_H_