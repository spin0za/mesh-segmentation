#ifndef  _DYNAMIC_DOUBLE_COVERING_MESH_H_
#define  _DYNAMIC_DOUBLE_COVERING_MESH_H_

#include <map>
#include <vector>
#include <queue>

#include "Mesh/Vertex.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/boundary.h"
#include "Mesh/iterators.h"
#include "Mesh/BaseMesh.h"
#include "Mesh/DynamicMesh.h"


namespace MeshLib
{

namespace Topology
{
/*-------------------------------------------------------------------------------------------

 Vertex Trait 

--------------------------------------------------------------------------------------------*/
/*! \brief CDynamicDoubleCoveringVertex class
 *
*    Vertex class for dynamic double covering
*/
class CDynamicDoubleCoveringVertex : public  CVertex
{
public:
	/*! CDynamicDoubleCoveringVertex constructor */
  CDynamicDoubleCoveringVertex() { m_father = 0; m_boundary = false; m_touched = false; };
  /*! CDynamicDoubleCoveringVertex destructor */
  ~CDynamicDoubleCoveringVertex() {};
  /*! father */
  int & father() { return m_father; };
  /*! dual */
  CDynamicDoubleCoveringVertex *& dual() { return m_dual; };
  bool & touched() { return m_touched; };
  CPoint2 & huv()  { return m_uv; };
	
  /*!
   *	vertex color
   */
  CPoint & rgb()  { return m_rgb; };

  void _to_string();

protected:
	/*! father trait */
	int m_father;
	/*! touched */
	bool m_touched;
	/*! uv */
	CPoint2 m_uv;
	/*! rgb */
	CPoint  m_rgb;
	/*! dual */
	CDynamicDoubleCoveringVertex * m_dual;


};

inline void CDynamicDoubleCoveringVertex::_to_string()
{
	CParser parser( m_string );
	parser._removeToken( "rgb" );
	
	parser._toString( m_string );
	
	std::stringstream iss;
	
	iss << "rgb=(" << m_rgb[0] << " " << m_rgb[1] << " " << m_rgb[2] << ")";

	if( m_string.size() > 0 )
	{
	  m_string += " ";
	}
	m_string += iss.str();
};

/*-------------------------------------------------------------------------------------------

 Edge Trait 

--------------------------------------------------------------------------------------------*/
/*! \brief CDynamicDoubleCoveringEdge class
 *
*    Edge class for dynamic double covering
*/
class CDynamicDoubleCoveringEdge : public  CEdge
{
public:
	/*! CDynamicDoubleCoveringEdge constructor */
  CDynamicDoubleCoveringEdge() { m_father = 0; m_sharp = false; };
  /*! CDynamicDoubleCoveringEdge destructor */
  ~CDynamicDoubleCoveringEdge() {};
  /*! father */
  int & father() { return m_father; };
  /*! boundary */
  bool & sharp() {return m_sharp; };
  /*! edge length */
  double & length() { return m_length; };
  /*! edge index */
  int  & index()  { return m_index; };
  /*  output to string */
  void _to_string();

protected:
	/*! father trait */
	int m_father;
	/*! boundary */
	bool m_sharp;
	/*! edge length */
	double m_length;
	/*! index */
	int    m_index;
};

inline void CDynamicDoubleCoveringEdge::_to_string()
{
		CParser parser( m_string );
		parser._removeToken( "l" );
		parser._toString( m_string );
		

		std::stringstream iss;
		iss << "l=("<< std::setprecision(12) << m_length << ")";
		if( m_string.length() > 0 )
		{
			m_string += " ";
		}
		m_string += iss.str();

}

/*-------------------------------------------------------------------------------------------

 Face Trait 

--------------------------------------------------------------------------------------------*/
/*! \brief CDynamicDoubleCoveringEdge class
 *
*    Edge class for dynamic double covering
*/
class CDynamicDoubleCoveringFace : public  CFace
{
public:
	/*! CDynamicDoubleCoveringFace constructor */
  CDynamicDoubleCoveringFace() { m_father = 0; };
  /*! CDynamicDoubleCoveringFace destructor */
  ~CDynamicDoubleCoveringFace() {};
  /*! father */
  int & father() { return m_father; };
  /*! touched */
  bool & touched() { return m_touched; };

protected:
	/*! father trait */
	int m_father;
	/*! touched */
	bool m_touched;
};

/*-------------------------------------------------------------------------------------------

 Half Edge Trait 

--------------------------------------------------------------------------------------------*/
/*! \brief CDynamicDoubleCoveringHalfEdge class
 *
*    HalfEdge class for dynamic double covering
*/
class CDynamicDoubleCoveringHalfEdge : public  CHalfEdge
{
public:
	/*! CDynamicDoubleCoveringHalfEdge constructor */
  CDynamicDoubleCoveringHalfEdge() { m_father = 0; };
  /*! CDynamicDoubleCoveringHalfEdge destructor */
  ~CDynamicDoubleCoveringHalfEdge() {};
  /*! father */
  int & father() { return m_father; };
protected:
	/*! father trait */
	int m_father;
};




/*! \brief CDoubleCoveringMesh class
*   
*   Compute the double covering of a mesh with boundaries
*   preassumptions
*	1. no edge is a loop, namely, the two end vertices of any edge are different
*   2. for any face, at least one vertex is not on the boundary
*
*/
template<typename V, typename E, typename F, typename H>
class CDynamicDoubleCoveringMesh : public CDynamicMesh<V,E,F,H>
{
public:
	typedef V CVertex;
	typedef E CEdge;
	typedef F CFace;
	typedef H CHalfEdge;

	typedef MeshVertexIterator<V,E,F,H>        MeshVertexIterator;
	typedef MeshFaceIterator<V,E,F,H>          MeshFaceIterator;
	typedef MeshEdgeIterator<V,E,F,H>          MeshEdgeIterator;
	typedef FaceVertexIterator<V,E,F,H>        FaceVertexIterator;
	typedef FaceHalfedgeIterator<V,E,F,H>	   FaceHalfedgeIterator;
	/*!
	 *	Compute the double covering
	 */
	void DoubleCovering();

};






template<typename V, typename E, typename F, typename H>
void CDynamicDoubleCoveringMesh<V,E,F,H>::DoubleCovering()
{
	write_vef( "temp.vef" );

	//--record the original vertices
	std::list<V *>origvertices;
	int vmax=0;
	for(MeshVertexIterator viter(this); !viter.end(); viter++)
	{
		V * v=*viter;
		origvertices.push_back(v);
		vmax = (vmax < v->id())? v->id():vmax;

		std::string line;
		std::stringstream iss(line);
		iss << "father=("<<v->id()<<")";
		
		if( v->string().length() > 0 )
		{
			v->string() += " ";
		}
		v->string() += iss.str();
	}

	//--record the original edges
	std::map<int, E*> m_edge_map;
	std::list<E*>origedges;
	int emax=0;
	for(MeshEdgeIterator eiter(this); !eiter.end(); eiter++)
	{
		E * e=*eiter;
		origedges.push_back(e);
		emax = (emax < e->id())? e->id():emax;

		std::string line;
		std::stringstream iss(line);
		iss << "father=("<<e->id()<<")";
		
		if( e->string().length() > 0 )
		{
			e->string() += " ";
		}
		e->string() += iss.str();

		m_edge_map.insert( std::pair<int,E*>(e->id(),e) );
	}


	//--record the original faces
	std::list<F *>origfaces;
	int fmax=0;
	for(MeshFaceIterator fiter(this); !fiter.end(); fiter++)
	{
		F * f=*fiter;

		std::string line;
		std::stringstream iss(line);
		iss << "father=("<<f->id()<<")";

		if( f->string().length() > 0 )
		{
			f->string() += " ";
		}
		f->string() += iss.str();

		origfaces.push_back(f);
		fmax = ( fmax < f->id() )? f->id():fmax;

	}


	std::map<int, int>corres;
	int currentid=vmax+1;
	for(std::list<V *>::iterator viter=origvertices.begin(); viter!=origvertices.end(); viter++)
	{
		V *v=*viter;
		if(!isBoundary(v))
		{
			V *newv=createVertex(currentid);
			newv->point()=v->point();
			newv->string() = v->string();
			corres[v->id()]=currentid;
			++currentid;
/*
			std::string line;
			std::stringstream iss(line);
			iss << "father=("<<v->id()<<")";
			
			if( newv->string().length() > 0 )
			{
				newv->string() += " ";
			}
			newv->string() += iss.str();
*/
		}
		else
			corres[v->id()]=v->id();
	}


	std::map<int, int>ecorres;
	int currenteid=emax+1;
	for(std::list<E *>::iterator eiter=origedges.begin(); eiter!=origedges.end(); eiter++)
	{
		E *e=*eiter;
		if(!isBoundary(e))
		{
			E *newe= new E;
			m_edges.push_back( newe );
			assert( newe != NULL );
			newe->id() = currenteid;	
			newe->string() = e->string();
			ecorres[e->id()]=currenteid;
			m_edge_map.insert( std::pair<int,E*>(newe->id(),newe ) );
			++currenteid;
/*
			std::string line;
			std::stringstream iss(line);
			iss << "father=("<<e->id()<<")";
			
			if( newe->string().length() > 0 )
			{
				newe->string() += " ";
			}
			newe->string() += iss.str();
*/
			H * ph[2];
			H * oldph[2];

			for( int k = 0; k < 2; k ++ )
			{
				//generate a new halfedge
				oldph[k] = edgeHalfedge(e,k);
				ph[k] = new H;
				assert( ph[k] );
				//link the new halfedge to the new edge
				newe->halfedge(k) = ph[k];
				ph[k]->edge() = newe;
				//link the new halfedge to the target vertex
				V * oldv = halfedgeVertex(oldph[k]);
				V * newv = idVertex( corres[vertexId( oldv )]);
				ph[k]->vertex() = newv;
				newv->halfedge() = ph[k];

				//the face link is null
				ph[k]->face() = NULL;
			}

		}
		else
		{
			ecorres[e->id()]=e->id();
			H * oldph = edgeHalfedge( e,0 );
			H * ph = new H;
			//link the new halfedge to the boundary edge
			e->halfedge(1) = ph;
			ph->edge() = e;
			//link the new halfedge to the target vertex
			V * pv = halfedgeSource( oldph );
			ph->vertex() = pv;
			pv->halfedge() = ph;
			//the face link is null
			ph->face() = NULL;
		}
	}


	for(std::list<F *>::iterator fiter=origfaces.begin(); fiter!=origfaces.end(); fiter++)
	{
		F * f=*fiter;

		H * oldh[3];
		H * newh[3];
		E * olde[3];
		E * newe[3];

		int i=0;
		for(FaceHalfedgeIterator hiter(f); !hiter.end(); hiter++)
		{
			oldh[i] =*hiter;
			olde[i] = halfedgeEdge( oldh[i] );
			int neweid = ecorres[olde[i]->id()];
			newe[i] = m_edge_map[neweid];
			i ++;
		}
		
		for( int i = 0; i < 3; i ++ )
		{
			V * oldv = halfedgeTarget( oldh[i] );
			V * newv = idVertex( corres[oldv->id()] );

			H * ph = edgeHalfedge( newe[i],0);
			if( ph->vertex() == newv )
			{
				ph = edgeHalfedge( newe[i],1);
			}
			newh[i]  = ph;
		}

		//convert orientation
		H * th = newh[0];
		newh[0] = newh[1];
		newh[1] = th;

		for( int i = 0; i < 3; i ++ )
		{
			newh[i]->he_next()  = newh[(i+1)%3];
			newh[i]->he_prev()  = newh[(i+2)%3];
		}

		F * newf = new F;
		assert( newf != NULL );
		m_faces.push_back( newf );
		newf->id() = fmax + 1;
		fmax ++;
		newf->string() = f->string();

		newf->halfedge() = newh[0];
		for( int i = 0; i < 3; i ++ )
		{
			newh[i]->face()  = newf;
		}

		m_map_face.insert( std::pair<int,F*>( newf->id(), newf ) );

	}

	for( std::list<E*>::iterator eiter = m_edges.begin(); eiter != m_edges.end(); eiter ++ )
	{
		E * e = *eiter;
		H * ph[2];
		for( int i = 0; i < 2; i ++ ) ph[i] = edgeHalfedge( e,i );
		if( ph[0]->vertex()->id() < ph[1]->vertex()->id() )
		{
			e->halfedge(0) = ph[1];
			e->halfedge(1) = ph[0];
		}
	}
	

};



typedef CDynamicDoubleCoveringMesh<CDynamicDoubleCoveringVertex, CDynamicDoubleCoveringEdge, CDynamicDoubleCoveringFace, CDynamicDoubleCoveringHalfEdge> CDDCMesh;

} //namespace Topology
} //namespace MeshLib
#endif