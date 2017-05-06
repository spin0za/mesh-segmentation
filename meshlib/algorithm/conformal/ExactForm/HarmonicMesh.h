/*!
*      \file HarmonicMesh.h
*      \brief Mesh for Computing Exact Harmonic Forms
*	   \author David Gu
*      \date Document 10/11/2010
*
*/

#ifndef  _HARMONIC_MESH_H_
#define  _HARMONIC_MESH_H_

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
#include "Parser/traits_io.h"

namespace MeshLib
{

namespace Holomorphy
{
/*!
 * \brief CHVertex class
 *
 * vertex class for computing harmonic exact forms
 * adding vertex index trait, vertex function u trait, and vertex uv field
 */

class CHVertex : public CVertex
{
public:
	/*! CHVertex constructor
	*/
	CHVertex() { m_index = 0;  };
	/*! CHVertex destructor
	 */
	~CHVertex() {};
	/*! Vertex index
	*/
	int &    idx() { return m_index; };
	/*! Vertex texture coordinates
	*/
	CPoint2 & uv() { return m_uv; };
	/*!	Vertex function
	*/
	double  &  u() { return m_u;  };
	/*!	Vertex z
	 */
	std::complex<double> & z() { return m_z; };

	/*!	Vertex mu
	 */
	std::complex<double> & mu() { return m_mu; };
	/*!	
	 * whether the vertex is fixed or not
	 */
	bool & fixed() { return m_fixed; };

	/*! Convert vertex texture coordinates to string
	 */
	void  _to_string()
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

protected:

	/*! Vertex index */
	int         m_index;
	/*! Vertex texture coordinates */
	CPoint2     m_uv;
	/*! Vertex function */
    double      m_u;
	/*! Vertex mu */
	std::complex<double> m_mu;
	/*! Vertex z */
	std::complex<double> m_z;
	/*! fixed */
	bool m_fixed;
};

/*! \brief CHEdge class
 *
 *  Edge class for computing exact harmonic form,
 *  add edge weight trait, and 1-form on edge trait.
 */

class CHEdge : public  CEdge
  {
  public:
    /*! CHEdge constructor
	 */
	 CHEdge() { m_weight = 0; m_du = 0; };
    /*! CHEdge destructor
	 */
    ~CHEdge(){};

	/*!	Edge 1-form
	 */
	double & du() { return m_du; };
	/*! Edge weight
	 */
	double & weight() { return m_weight; };
	/*! Edge length
	 */
	double & length() { return m_length; };
	
	/*! Convert edge 1-form du to string.
	*/
	void _to_string()
	{
		CParser parser( m_string );
		parser._removeToken( "du" );

		parser._toString( m_string );
		
		std::string line;
		std::stringstream iss(line);
		iss << "du=(" << m_du << ")";

		if( m_string.length() > 0 )
		{
			m_string += " ";
			m_string += iss.str();
		}
		else
		{
			m_string = iss.str();
		}
	};

  protected:
	  /*! Edge weight */
	double   m_weight;
	/*! Edge 1-form */
    double   m_du;
	/*! Edge length */
	double   m_length;
};
	
/*! \brief CHHalfEdge class
* 
*	Halfedge Class for computing exact harmonic form
*/
class CHHalfEdge : public  CHalfEdge
  {

  public:
	/*! CHHalfEdge constructor
	*/
	 CHHalfEdge() { m_angle = 0; };
	/*! CHHalfEdge destructor
	*/
    ~CHHalfEdge(){};
	  /*! Corner angle */
	double &   angle() { return m_angle; };
	
	//void _from_string();
	//void _to_string();

  protected:
	  /*! Corner angle */
	  double   m_angle;

};

/*! \brief CHMesh class
 *
 *   Mesh class for computing exact harmonic form
 */
template<typename V, typename E, typename F, typename H>
class CHMesh : public CBaseMesh<V,E,F,H>
{
public:
	typedef V CVertex;
	typedef E CEdge;
	typedef F CFace;
	typedef H CHalfEdge;

	typedef CBoundary<V,E,F,H>  CBoundary;
	typedef CLoop<V,E,F,H> CLoop;
	typedef MeshVertexIterator<V,E,F,H> MeshVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
	typedef MeshEdgeIterator<V,E,F,H>   MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef MeshFaceIterator<V,E,F,H>     MeshFaceIterator;
	
public:
};


typedef CHMesh<CHVertex, CHEdge, CFace, CHHalfEdge> CHarmonicMesh;

} //namespace Holomorphy

unsigned long long Holomorphy::CHarmonicMesh::m_input_traits  = 0;
unsigned long long Holomorphy::CHarmonicMesh::m_output_traits = VERTEX_UV;

} //namespace Holomorphy
#endif  _HARMONIC_MESH_H_