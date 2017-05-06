/*!
*      \file ZipperMesh.h
*      \brief Mesh for computer conformal mapping using zipper algorithm
*	   \author David Gu
*      \date Documented 02/13/2014
*
*/

#ifndef  _ZIPPER_MESH_H_
#define  _ZIPPER_MESH_H_

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
/*!
*	\brief CZipperVertex class
*
*	Vertex class for Zipper
*/
class CZipperVertex : public CVertex
{

public:
	/*!
	 *	CHarmonicVertex constructor
	 */
	CZipperVertex() {m_k=0; };
	/*!
	 *	CHarmonicVertex destructor
	 */
	~CZipperVertex() {};
	/*
	 *	Zipper curvature
	 */
	double & k() { return m_k; };
	/*
	 *	Vertex normal
	 */
	CPoint & normal() { return m_normal; };
	/*
	 *	Complex coordinates
	 */
	std::complex<double> & z() { return m_z; };
	/*
	 *	to string
	 */
	void _to_string();

protected:
	/*
	 *	Zipper curvature
	 */
	double m_k;
	/*!
	 *	Vertex normal
	 */
	CPoint m_point;
	/*!
	 *	Complex coordinates
	 */
	std::complex<double> m_z;
};

/*! write v->huv to the string 
*/
inline void CZipperVertex::_to_string()
{
	CParser parser( m_string );
	parser._removeToken( "uv" );
	parser._toString( m_string );
	std::stringstream iss;
	iss << "uv=(" << m_z.real() << " " << m_z.imag() << ")";

	if( m_string.length() > 0 )
	{
		m_string += " ";
	}
	m_string += iss.str();
}

/*!
*	\brief CZipperEdge class
*
*	Edge class for Zipper curvature
*/
class CZipperEdge : public  CEdge
  {
  public:
    /*!	CZipperEdge constructor
	 */
	 CZipperEdge() { m_length=0;  };
    /*!	CZipperEdge destructor
	 */
    ~CZipperEdge(){};
	/*! edge length trait
	 */
	double & length() { return m_length; };

  protected:
	/*! edge length trait */
	double   m_length;
};



/*!
*	\brief CZipperEdge class
*
*	Face class for Zipper curvature
*/
class CZipperFace : public  CFace
  {
  public:
    /*!	CZipperFace constructor
	 */
	 CZipperFace() {};
    /*!	CZipperFace destructor
	 */
    ~CZipperFace(){};
	/*!
	 *	face area
	 */
	double & area() { return m_area; };
	/*!
	 *	face normal
	 */
	CPoint & normal(){ return m_normal; };
 
  protected:
	double m_area;
	CPoint m_normal;
};


/*!
*	\brief CZipperHalfEdge class
*
*	HalfEdge class for Zipper curvature
*/
class CZipperHalfEdge : public  CHalfEdge
  {
  public:
    /*!	CHarmonicHalfEdge constructor
	 */
	 CZipperHalfEdge() {};
    /*!	CHarmonicHalfEdge destructor
	 */
    ~CZipperHalfEdge(){};
	/*!	Corner angle trait
	 */
	double & angle() { return m_angle; };

  protected:
	  /*! Corner angle trait */
	double m_angle;
};

/*-------------------------------------------------------------------------------------------------------------------------------------

	ZipperCurvature Mesh Class

--------------------------------------------------------------------------------------------------------------------------------------*/
/*!
 *	\brief CZipperMesh class
 *
 *	Mesh class for Zipper curvature
 */
template<typename V, typename E, typename F, typename H>
class CZipperMesh : public CBaseMesh<V,E,F,H>
{
public:
	
	typedef V CVertex;
	typedef E CEdge;
	typedef F CFace;
	typedef H CHalfEdge;

	typedef CZipperMesh<V,E,F,H> M;

	typedef CBoundary<V,E,F,H> CBoundary;
	typedef CLoop<V,E,F,H> CLoop;
	
	typedef MeshVertexIterator<V,E,F,H> MeshVertexIterator;
	typedef MeshEdgeIterator<V,E,F,H> MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
	typedef MeshFaceIterator<V,E,F,H> MeshFaceIterator;
	typedef FaceVertexIterator<V,E,F,H> FaceVertexIterator;
	typedef VertexFaceIterator<V,E,F,H> VertexFaceIterator;
	typedef FaceHalfedgeIterator<V,E,F,H> FaceHalfedgeIterator;
	typedef VertexInHalfedgeIterator<V,E,F,H> VertexInHalfedgeIterator;
};

/*! Mesh class for CZipper class, Abbreviated as 'CZMesh'
 */
typedef CZipperMesh<CZipperVertex, CZipperEdge, CZipperFace, CZipperHalfEdge> CZMesh;	

};
#endif  _ZIPPER_MESH_H_