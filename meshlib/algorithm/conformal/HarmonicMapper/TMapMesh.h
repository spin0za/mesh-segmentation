/*!
*      \file TMapMesh.h
*      \brief Mesh for Computing TMap 
*	   \author David Gu
*      \date Document 09/29/2013
*
*/
#ifndef  _TMAP_MESH_H_
#define  _TMAP_MESH_H_

#include <map>
#include <vector>
#include <queue>
#include <complex>

#include "Mesh/BaseMesh.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "mesh/boundary.h"
#include "Parser/parser.h"
#include "Parser/traits_io.h"

namespace MeshLib
{

namespace Holomorphy
{

/*!
*	\brief CELVertex class
*
*	Vertex class for extremal length
*/
class CTMapVertex : public CVertex
{
public:

	/*!
	 *	CELVertex constructor
	 */
	CTMapVertex() { m_index = 0; m_valence = 0; };
	/*!
	 *	CELVertex destructor
	 */
	~CTMapVertex() {};

	/*! 
	 *	Vertex index
	 */
	int &     idx()    { return m_index;  };
	
	/*! Topological valence
	 */
	int & valence()    { return m_valence; };

	/*! whether the vertex is fixed
	 */
	bool & fixed()    { return m_fixed; };

	/*!	Vertex texture uv coordinates
	 */
	CPoint2 & uv() { return m_uv; };
	/*!	Vertex function u
	*/
	double  &  u() { return m_u;  };
	/*! convert vertex traits to the vertex string */
	void  _to_string();
	/*!
		conformal parameter
	*/
	std::complex<double> & z() { return m_z; };
	/*!
		image complex parameter
	*/
	std::complex<double> & w() { return m_w; };


protected:	//output
	/*! vertex texture coordinates */
	CPoint2     m_uv;	
protected:
	/*! vertex index */
	int         m_index;
	/*! vertex topological valence */
	int         m_valence;
	/*! whether the current vertex is fiexed */
	bool        m_fixed; 
	/*! vertex function */
    double      m_u;
	/*! complex parameter */
	std::complex<double> m_z;
	/*! complex parameter for the image*/
	std::complex<double> m_w;
};


/*! write vertex traits uv to the string */
inline	void  CTMapVertex::_to_string()
	{
		CParser parser( m_string );
		parser._removeToken( "uv" );

		parser._toString( m_string );
		
		std::string line;
		std::stringstream iss(line);
		iss << "uv=(" << m_uv[0] << " " << m_uv[1] << ")";

		if( m_string.length() > 0 )
		{
			m_string += " ";
		}
		m_string += iss.str();

	};
/*! \brief CTMapEdge class
*
* Edge class for computing extremal length
*/

class CTMapEdge : public  CEdge
  {
  public:
  /*! CTMapEdge constructor
  */
	 CTMapEdge() { m_weight = 0; m_sharp = false; };
  /*! CELEdge destructor
  */
    ~CTMapEdge(){};
	/*! Edge weight */
	double & weight() { return m_weight; };
	/*! Edge length */
	double & length() { return m_length; };


	/*! sharp flag for the edge */
	bool & sharp() { return m_sharp; };

  protected: //output
	/*! edge weight */
	double   m_weight;
	/*! edge sharp flat */
	bool     m_sharp;
	/*! edge length */
	double   m_length;
};



/*! \brief CTMapHalfEdge class
*
*  Halfedge class for computing extremal length
*/
class CTMapHalfEdge : public  CHalfEdge
  {
  public:
	/*! CHCFHalfEdge constructor */
	 CTMapHalfEdge() { m_angle = 0; };
	/*! CHCFHalfEdge destructor */
    ~CTMapHalfEdge(){};
	  /*! corner angle */
	double &   angle() { return m_angle; };
	
  //no output
  protected: 
	  /*! corner angle */
	  double   m_angle;

};


/*! \brief CTapFace class
*
*  Halfedge class for computing extremal length
*/
class CTMapFace : public  CFace
  {
  public:
	/*! CTMapFace constructor */
	 CTMapFace() {};
	/*! CTMapFace destructor */
    ~CTMapFace(){};
	  /*! face mu */
	std::complex<double> &   mu() { return m_mu; };
	  /*! face nu */
	std::complex<double> &   nu() { return m_nu; };
	
  //no output
  protected: 
	  /*! face Beltrami coefficient */
	  std::complex<double>   m_mu;
	  /*! face Beltrami coefficient */
	  std::complex<double>   m_nu;

};
/*! \brief CTMapMesh class
*
*  Mesh class for computing TMap 
*/
template<typename V, typename E, typename F, typename H>
class CTMapMesh : public CBaseMesh<V,E,F,H>
{
public:
	typedef V CVertex;
	typedef E CEdge;
	typedef F CFace;
	typedef H CHalfEdge;

	typedef CBoundary<V,E,F,H>  CBoundary;
	typedef CLoop<V,E,F,H> CLoop;
	typedef CLoopSegment<V,E,F,H> CLoopSegment;

	typedef MeshVertexIterator<V,E,F,H> MeshVertexIterator;
	typedef MeshFaceIterator<V,E,F,H> MeshFaceIterator;
	typedef MeshEdgeIterator<V,E,F,H> MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
	typedef VertexOutHalfedgeIterator<V,E,F,H> VertexOutHalfedgeIterator;
	typedef FaceHalfedgeIterator<V,E,F,H> FaceHalfedgeIterator;
	typedef FaceVertexIterator<V,E,F,H> FaceVertexIterator;
public:
};

typedef CTMapMesh<CTMapVertex, CTMapEdge, CTMapFace, CTMapHalfEdge> CTMAPMesh;

} //namespace Holomorphy
} //namespace MeshLib

#endif  _TMAP_MESH_H_