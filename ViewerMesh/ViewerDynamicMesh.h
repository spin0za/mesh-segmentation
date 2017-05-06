/*! \file ViewerDynamicMesh.h
*   \brief Viewer DynamicMesh
*   \author David Gu
*   \date   documented on 10/07/2013
*
*   Dynamic Mesh for viewer 
*/
#ifndef  _VIEWER_DYNAMIC_MESH_H_
#define  _VIEWER_DYNAMIC_MESH_H_

//#define INFINITY 1e6

#include <map>
#include <vector>

#include "Mesh/Vertex.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "Mesh/HalfEdge.h"

#include "Mesh/DynamicMesh.h"
#include "Mesh/boundary.h"
#include "Mesh/iterators.h"
#include "Parser/parser.h"
#include "Parser/traits_io.h"

namespace MeshLib
{
/*! \brief CViewerHalfEdge class
*
*	HalfEdge class for viewer
*   Trait: child
*/
  class CViewerHalfEdge: public CHalfEdge
  {
  protected:
    /*! child */
    CViewerHalfEdge * m_child;
	/* angle at target vertex */
	double m_angle;

  public:
    /*! face normal */
	  CViewerHalfEdge* & child() { return m_child; };
	  /* angle at target vertex */
	  double & angle() { return m_angle; };
  };

/*! \brief CViewerFace class
*
*	Face class for viewer
*   Trait: face normal
*/
  class CViewerFace: public CFace
  {
  protected:
	/*! face color */
	CPoint m_rgb;
    /*! face normal */
    CPoint m_normal;
	CPoint m_normal_filtered;
	/*! if face is picked*/
	bool   m_picked;
	/* area */
	double m_area;
	/* obtuse vertex*/
	int m_obtuseVertex;
	/* metric used in region growing */
	double m_metric;
	/* mark */
	int  m_patch_label;
	bool m_isCandidate;
	bool m_boundary_mark;
	int	 m_isomark;
	int  m_layer;
	CViewerFace * m_parent;
	CViewerFace * m_marked_neighbor;
  public:
	  CViewerFace()
	  {
		  m_rgb = CPoint(1, 1, 1); //default color is white
		  m_normal_filtered = CPoint(0, 0, 0);
		  m_picked = false;
		  m_area = 0;
		  m_obtuseVertex = -1;
		  m_patch_label = -1;
		  m_isCandidate = false;
		  m_boundary_mark = false;
		  m_isomark = -1;
		  m_layer = 0;
		  /* metric used in region growing */
		  m_metric = INFINITY;
		  m_parent = nullptr;
		  m_marked_neighbor = nullptr;
	  }
	  /*! face normal */
	  CPoint & normal() { return m_normal; };
	  CPoint & normal_filtered() { return m_normal_filtered; };
	  /*! face color */
	  CPoint  & rgb()	{ return m_rgb; };
	  /*! picked */
	  bool   & picked() { return m_picked; };
	  /* area */
	  double & area() { return m_area; };
	  /* obtuse vertex*/
	  int & obtuseVertex() { return m_obtuseVertex; };
	  /* mark */
	  int & patch_label() { return m_patch_label; };
	  bool & isCandidate() { return m_isCandidate; };
	  bool & boundary_mark() { return m_boundary_mark; };
	  int & isomark() { return m_isomark; };
	  int & layer() { return m_layer; };
	  CViewerFace *& parent() { return m_parent; };
	  /* metric used in region growing */
	  double & metric() { return m_metric; };
	  bool operator < (const CViewerFace *&f) { return m_metric >= f->m_metric; };
	  CViewerFace *& marked_neighbor() { return m_marked_neighbor; };
  };

  /*! \brief CViewerVertex class
  *   
  *   Vertex class for viewer
  *   Trait : vertex rgb color
  */
  class CViewerVertex : public  CVertex
  {
  protected:
	  /*! vertex color */
	  CPoint	 m_yuv;
	  CPoint	 m_rgb;
	  CPoint	 m_rgb_old;
	  CPoint	 m_rgb_bak;
	  /*! normal map */
	  CPoint   m_normal_map;
	  /*! father */
	  int      m_father;
	  /*! is current vertex picked */
	  bool     m_picked;
	  /* scalar field */
	  double   m_field;
	  /* mark */
	  int	 m_part_label;
	  int	 m_part_label_prev;
	  int	 m_isomark;
	  bool m_colormark;
	  /* mean curvature*/
	  double   m_meanCurvature;
	  /* standard */
	  int m_sid;
	  int m_degree;
  public:
	  /*! CViewerVertex Constructor */
	  CViewerVertex()
	  {
		  m_yuv = CPoint(1, 0.5, 0.5);
		  m_rgb = CPoint(1, 1, 1); //default color is white
		  m_rgb_old = CPoint(1, 1, 1);
		  m_rgb_bak = CPoint(1, 1, 1);
		  m_normal_map = CPoint(1,0,0);
		  m_picked = false;

		  m_part_label = -1;
		  m_part_label_prev = -1;
		  m_isomark = -1;
		  m_colormark = false;
		  m_field = 1;
		  m_meanCurvature = 0;
		  m_sid = -1;
		  m_degree = 0;
	  }
	  /*! vertex rgb color */
	  CPoint  & yuv()	{ return m_yuv;	};
	  CPoint  & rgb()	{ return m_rgb; };
	  CPoint  & rgb_old()	{ return m_rgb_old; };
	  CPoint  & rgb_bak()	{ return m_rgb_bak; };
	
	  /*! normal map */
	  CPoint  & normal_map() { return m_normal_map; };

	  /*! father */
	  int & father() { return m_father; };

	  /*! read vertex rgb, uv from vertex string */
	  void _from_string();

	  /*! if the current vertex is selected */
	  bool & picked() { return m_picked; };

	  /* mark */
	  int & part_label() { return m_part_label; };
	  int & part_label_prev() { return m_part_label_prev; };
	  int & isomark() { return m_isomark; };
	  bool & colormark() { return m_colormark; };

	  /* scalar field */
	  double & field() { return m_field; };

	  /* mean curvature */
	  double & meanCurvature() { return m_meanCurvature; };
	  /* standard */
	  int & sid() { return m_sid; };
	  int & degree() { return m_degree; };
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
			  else
			  if( token->m_key == "normal" )
			  {
				  token->m_value >> m_normal_map;
			  }
			  else
			  if( token->m_key == "father" )
			  {
				  std::string line = strutil::trim( token->m_value, "()");
				  m_father = strutil::parseString<int>( line );	
			  }
		  }
  };


  /*! \brief CViewerEdge class
  *   
  *   Edge class for viewer
  *   Trait : Edge sharp
  */
  class CViewerEdge : public  CEdge
  {
  protected:
	  /*! edge sharp */
	  bool   m_sharp;
	  /* mark */
	  int	 m_isomark;
	  /* intersection with isoline */
	  CPoint m_iso_cross;
	  CPoint m_bd_cross;
	  /* standard id */
	  int m_sid;
	  /* similarity */
	  double m_similarity;
	  /* weight */
	  double m_weight[2];
	
  public:
	  /*! CViewerVertex Constructor */
	  CViewerEdge()
	  {
		  m_sharp = false;
		  m_isomark = -1;
		  m_iso_cross = CPoint(0, 0, 0);
		  m_bd_cross = CPoint(0, 0, 0);
		  m_sid = -1;
		  m_similarity = -1;
		  m_weight[0] = -1;
		  m_weight[1] = -1;
	  }
	  /*! vertex rgb color */
	  bool  & sharp()    { return m_sharp; };
	  /* mark */
	  int & isomark()	{ return m_isomark; };
	  /* intersection with isoline */
	  CPoint & iso_cross() { return m_iso_cross; };
	  CPoint & bd_cross() { return m_bd_cross; };
	  /* standard id */
	  int & sid() { return m_sid; };
	  /* similarity */
	  double & similarity() { return m_similarity; };
	  /* weight */
	  double & weight(CViewerVertex *v)	  { return v == m_halfedge[0]->source() ? m_weight[0] : m_weight[1]; };

	  /*! read vertex rgb, uv from vertex string */
	  void _from_string();
  };
 
 // read vertex rgb, uv from vertex string 
 inline void CViewerEdge::_from_string()
 {
		  CParser parser( m_string );
		
		  for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
		  {
			  CToken * token = *iter;
			  if( token->m_key == "sharp" )
			  {
				  m_sharp = true;
			  }
			  break;
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
class CViewerDynamicMesh : public CDynamicMesh<V,E,F,H>
{
public:
	
	typedef V CVertex;
	typedef E CEdge;
	typedef F CFace;
	typedef H CHalfEdge;

	typedef CBoundary<V,E,F,H> CBoundary;
	typedef CLoop<V,E,F,H> CLoop;
	typedef MeshVertexIterator<V,E,F,H> MeshVertexIterator;
	typedef MeshFaceIterator<V,E,F,H>   MeshFaceIterator;
	typedef MeshEdgeIterator<V,E,F,H>   MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
	typedef VertexFaceIterator<V, E, F, H> VertexFaceIterator;
	typedef VertexInHalfedgeIterator<V, E, F, H> VertexInHalfedgeIterator;
	typedef VertexOutHalfedgeIterator<V, E, F, H> VertexOutHalfedgeIterator;
	typedef FaceHalfedgeIterator<V, E, F, H> FaceHalfedgeIterator;
	typedef FaceVertexIterator<V, E, F, H> FaceVertexIterator;

public:
	// the center of face f
	CPoint faceCenter(CViewerFace *f)
	{
		CViewerVertex *v[3];
		int vid = 0;
		for (CVDMesh::FaceVertexIterator viter(f); !viter.end(); viter++)
		{
			v[vid++] = *viter;
		}
		return (v[0]->point() + v[1]->point() + v[2]->point()) / 3;
	}

	// intersection of two edges
	CViewerVertex * eCommonVertex(CViewerEdge *e1, CViewerEdge *e2)
	{
		CViewerVertex * v11 = edgeVertex1(e1);
		CViewerVertex * v12 = edgeVertex2(e1);
		CViewerVertex * v21 = edgeVertex1(e2);
		CViewerVertex * v22 = edgeVertex2(e2);
		if (v11 == v21 || v11 == v22)
		{
			return v11;
		}
		else if (v12 == v21 || v12 == v22)
		{
			return v12;
		}
		else
		{
			return NULL;
		}
	};

	// the vertex in the edge e which is not v
	CViewerVertex * eOtherVertex(CViewerEdge *e, CViewerVertex *v)
	{
		if (v == edgeVertex1(e))
		{
			return edgeVertex2(e);
		}
		else if (v == edgeVertex2(e))
		{
			return edgeVertex1(e);
		}
		else
		{
			return NULL;
		}
	};

	// the face attached to the edge e which is not f
	CViewerFace * eOtherFace(CViewerEdge *e, CViewerFace *f)
	{
		if (f == edgeFace1(e))
		{
			return edgeFace2(e);
		}
		else if (f == edgeFace2(e))
		{
			return edgeFace1(e);
		}
		else
		{
			return NULL;
		}
	};

	// the vertex v in f which is opposite to e
	CViewerVertex * eOppoVertex(CViewerEdge *e, CViewerFace *f)
	{
		CViewerVertex *v1 = edgeVertex1(e);
		CViewerVertex *v2 = edgeVertex2(e);
		for (CVDMesh::FaceVertexIterator viter(f); !viter.end(); viter++)
		{
			CViewerVertex *v3 = *viter;
			if (v3 != v1 && v3 != v2)
			{
				return v3;
			}
		}
		return NULL;
	};

	// the vertex opposite to e in the face attached to e which is not f
	CViewerVertex * eNextVertex(CViewerEdge *e, CViewerFace *f)
	{
		CViewerFace *f1 = edgeFace1(e);
		CViewerFace *f2 = edgeFace2(e);
		if (f == f1)
		{
			return eOppoVertex(e, f2);
		}
		else
		{
			return eOppoVertex(e, f1);
		}
	};

	// the face defined by v and e
	CViewerFace * vEdgeFace(CViewerVertex *v, CViewerEdge *e)
	{
		CViewerFace *f1 = edgeFace1(e);
		if (vertexInFace(v, f1))
		{
			return f1;
		}
		CViewerFace *f2 = edgeFace2(e);
		if (vertexInFace(v, f2))
		{
			return f2;
		}
		return NULL;
	};

	// make sure e1's end is e2's head
	bool orientEdges(CViewerEdge *&e1, CViewerEdge *&e2)
	{
		CViewerHalfEdge *h10 = edgeHalfedge(e1, 0);
		CViewerHalfEdge *h11 = edgeHalfedge(e1, 1);
		CViewerHalfEdge *h20 = edgeHalfedge(e2, 0);
		CViewerHalfEdge *h21 = edgeHalfedge(e2, 1);
		if (halfedgeNext(h10) == h20 || halfedgeNext(h10) == h21 ||
			halfedgeNext(h11) == h20 || halfedgeNext(h11) == h21)
		{
			return true;
		}
		else if (halfedgeNext(h20) == h10 || halfedgeNext(h20) == h11 ||
				 halfedgeNext(h21) == h10 || halfedgeNext(h21) == h11)
		{
			CViewerEdge *e3 = e2;
			e2 = e1;
			e1 = e3;
			return true;
		}
		else
		{
			return false;
		}
	};

	bool vertexInFace(CViewerVertex *v0, CViewerFace *f)
	{
		for (CVDMesh::FaceVertexIterator viter(f); !viter.end(); viter++)
		{
			CViewerVertex *v = *viter;
			if (v == v0)
			{
				return true;
			}
		}
		return false;
	};

	// let the edge point from lower scalar to higher
	void tellDirection(CViewerEdge *e, CViewerVertex *&v1, CViewerVertex *&v2)
	{
		v1 = edgeVertex1(e);
		v2 = edgeVertex2(e);
		if (v1->field() > v2->field())
		{
			v1 = edgeVertex2(e);
			v2 = edgeVertex1(e);
		}
	};

	bool isPointOnEdge(CPoint p, CViewerEdge *e)
	{
		CViewerVertex *v1 = edgeVertex1(e);
		CViewerVertex *v2 = edgeVertex2(e);
		CPoint p1 = v1->point();
		CPoint p2 = v2->point();
		CPoint zero(0, 0, 0);
		CPoint cp = (p - p1) ^ (p2 - p);
		return (cp - zero).norm() < 1e-6;
	}
};

typedef CViewerDynamicMesh<CViewerVertex, CViewerEdge, CViewerFace, CViewerHalfEdge> CVDMesh;

unsigned long long CVDMesh::m_input_traits  = VERTEX_UV| VERTEX_RGB | EDGE_SHARP;
unsigned long long CVDMesh::m_output_traits = 0;

}
#endif  VIEWER_DYNAMIC_MESH_H_