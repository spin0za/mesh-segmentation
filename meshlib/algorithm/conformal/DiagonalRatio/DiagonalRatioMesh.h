/*!
*      \file DiagonalRatioMesh.h
*      \brief Mesh for computing holomorphic forms using diagonal ratio
*	   \author David Gu
*      \date Documented on 06/29/2011
*
*/
#ifndef  _DIAGONAL_RATIO_MESH_H_
#define  _DIAGONAL_RATIO_MESH_H_

#include <map>
#include <vector>
#include <complex>

#include "Mesh/BaseMesh.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "mesh/boundary.h"
#include "Parser/parser.h"

namespace MeshLib
{

namespace Holomorphy
{
 /*! \brief CDRVertex class
  *
  *  Vertex class for DiagonalRatioMesh
  *  Traits:
  *  vertex index : m_index
  *  vertex uv coordinates: m_uv
  */
  class CDRVertex : public  CVertex
  {
  protected:
	/*!  vertex index */
    int         m_index;
	/*!  vertex uv position*/
	CPoint2     m_huv;
	/*!  whether the vertex is fixed or not */
	bool		m_fixed;
  public:
	/*!  vertex index */
	int & idx()   { return m_index; };
	/*! vertex parameter*/
	CPoint2 & huv() { return m_huv; };
	/*! whether the vertex is fixed */
	bool & fixed()  { return m_fixed; };

  public:
	  /*! CDRVertex Constructor */
    CDRVertex() { m_index = 0; m_fixed = false; };
	/*! CDRVertex destructor */
    ~CDRVertex() {};
  };


/*--------------------------------------------------------------------------------------------------------------------------------------
	Diagonal Ratio Edge Trait

---------------------------------------------------------------------------------------------------------------------------------------*/
/*! \brief CDREdge class
* 
*	Edge class for Diagonal Ration Mesh
*   Traits:
*   Diagonal Ratio : rho, du
*/
class CDREdge : public  CEdge
  {
  protected:
    /*! Edge index */
    int  m_index;
	/*! Harmonic 1-form */
	double  m_du;
	/*! Diagonal Ratio */
	std::complex<double> m_rho;
	/*! Holomorphic 1-form */
	std::complex<double> m_duv;
	/*! Edge length */
	double m_length;
	/*! whether the edge is sharp or not */
	bool   m_sharp;

  public:
	  /*! CDREdge constructor */
    CDREdge() { m_sharp = false; };
	/*! CDREdge destructor */
    ~CDREdge(){};
	/*! read sharp trait from edge string */
	void _from_string();
	/*! write sharp trait from edge string */
	void _to_string();
	/*! 1-form of the edge */
	double & du()  { return m_du; };
	/*! Holomorphic 1-form */
	std::complex<double> & duv() { return m_duv; };
	/*! Diagonal Ratio */
	std::complex<double> & rho() { return m_rho; };
	/*! Edge index */
	int & idx() { return m_index; };
	/*! Edge length */
	double & length() { return m_length; };
	/*! Edge sharp  */
	bool   & sharp()  { return m_sharp;  };
};

/*! read sharp trait from edge string */

inline void CDREdge::_from_string()
{
	CParser parser( m_string );

	for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
	{
			CToken * token = *iter;
			if( token->m_key == "du" )
			{
				std::string line= strutil::trim( token->m_value, "()");
				m_du = strutil::parseString<double>( line );
				break;
			}
			if( token->m_key == "sharp" )
			{
				m_sharp = true;
			}
	};
};

/*! write sharp trait from edge string */

inline void CDREdge::_to_string()
{
	if( m_sharp )
	{
		if( m_string.length() > 0 )
		{
			m_string += " ";
		}
		m_string += "sharp";
	}
};


/*! \brief CDiagonalRatioMesh class
*
*	Mesh class for computing holomorphic 1-forms using diagonal ratio
*/
template<typename V, typename E, typename F, typename H>
class CDiagonalRatioMesh : public CBaseMesh<V,E,F,H>
{
public:
	typedef V							  CVertex;
	typedef E							  CEdge;
	typedef F							  CFace;
	typedef H							  CHalfEdge;

	typedef CBoundary<V,E,F,H>            CBoundary;
	typedef CLoop<V,E,F,H>                CLoop;
	typedef MeshVertexIterator<V,E,F,H>   MeshVertexIterator;
	typedef MeshFaceIterator<V,E,F,H>	  MeshFaceIterator;
	typedef MeshEdgeIterator<V,E,F,H>     MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H>	  VertexEdgeIterator;
	typedef FaceHalfedgeIterator<V,E,F,H> FaceHalfedgeIterator;
	typedef FaceVertexIterator<V,E,F,H>   FaceVertexIterator;
};


typedef CDiagonalRatioMesh<CDRVertex, CDREdge, CFace, CHalfEdge> CDRMesh;

} //namespace Holomorphy
} //namespace MeshLib

#endif  _DIAGONAL_RATIO_MESH_H_