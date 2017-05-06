/*! \file HolomorphicQuadraticFormMesh.h
*    \brief Mesh for Computing Holomorphic Quadratic Forms
*	 \author David Gu
*    \date Document 06/16/2013
*/

#ifndef  _HOLOMORPHIC_QUADRATIC_FORM_MESH_H_
#define  _HOLOMORPHIC_QUADRATIC_FORM_MESH_H_

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
/*-----------------------------------------------------------------------------------------------------------------------------------------

	Holomorphic Quadratic Form Face Trait


-------------------------------------------------------------------------------------------------------------------------------------------*/

/*!
 *	\brief CHoloQuadFormVertex class
 *
 *	Vertex class for computing holomorphic forms
 *  add vertex huv and mu
 */

class CHoloQuadFormVertex: public CVertex
{
public:
	/*! CHoloQuadFormVertex constructor */
	CHoloQuadFormVertex() {};
	/*! CHoloQuadFormVertex destructor  */
	~CHoloQuadFormVertex() {};
	
	/*! conformal parameter */
	std::complex<double> & z() { return m_z; };

protected:
	/*! conformal parameter */
	std::complex<double>  m_z;

public:
	/*! read harmonic 1-form from edge string */
	void _from_string();
};


//read holomorphic quadratic form vertex trait "uv" to the trait m_z
inline void CHoloQuadFormVertex::_from_string()
{
	CParser parser( m_string );
	for( std::list<CToken*>::iterator titer = parser.tokens().begin(); titer != parser.tokens().end(); titer ++ )
	{
		CToken * pT = *titer;
		if( pT->m_key == "uv" )
		{
			std::string line= strutil::trim( pT->m_value, "()");
			 CPoint2 uv;
			 pT->m_value >> uv;
			 m_z = std::complex<double>(uv[0],uv[1]);
			break;
		}
	}
};


/*!
 *	\brief CHoloQuadFormFace class
 *
 *	Face class for computing holomorphic quadratic forms
 */

class CHoloQuadFormFace: public CFace
{
protected:
	/*! face area trait */
	double  m_area;
	/*! normal to the face */
	CPoint  m_normal;

public:
	/*! CHoloFormFace constructor */
	CHoloQuadFormFace() { m_area = 0; };
	/*! CHoloFormFace destructor  */
	~CHoloQuadFormFace() {};

public:
	/*! Face area */
	double & area()   { return m_area;   };
	/*! normal vector to the face */
	CPoint & normal() { return m_normal; };
};


/*-----------------------------------------------------------------------------------------------------------------------------------------

	Holomorphic Form Edge Trait


-------------------------------------------------------------------------------------------------------------------------------------------*/

 /*! \brief CHoloQuadFormEdge class
  *
  *  Edge class for computing holomorphic quadratic forms
  *  add harmonic 1-form trait m_du
  *  holomorphic 1-form trait m_duv
  *  edge length m_length traits.
  */
class CHoloQuadFormEdge : public  CEdge
{
  protected:
	  /*! edge length */
	  double  m_length;
	  /*! harmonic 1-form */
	  double  m_du;
	  /*! holomorphic 1-form */
	  CPoint2 m_duv;
	  /*! holomorphic quadratic forms */
	  std::complex<double> m_dz;

  public:
	  /*! CHoloFormEdge constructor */
      CHoloQuadFormEdge() { m_length = 0.0; m_du = 0;  };
	  /*! CHoloFormEdge destructor */
      ~CHoloQuadFormEdge(){};
	  /*! length of the edge */
	  double  & length() { return m_length; };
	  /*! harmonic 1-form */
	  double  & du()     { return m_du;     };
	  /*! holomorphic 1-form */
	  CPoint2 & duv()    { return m_duv;    };
	  
	  /*! holomorphic quadratic form */
	  std::complex<double> & dz()    { return m_dz; };

	  /*! read harmonic 1-form from edge string */
	  void _from_string();
	  /*! write holomorphic 1-form to edge string */
	  void _to_string();
};

//read harmonic 1-form trait "du" to the trait m_du
inline void CHoloQuadFormEdge::_from_string()
{
	CParser parser( m_string );
	for( std::list<CToken*>::iterator titer = parser.tokens().begin(); titer != parser.tokens().end(); titer ++ )
	{
		CToken * pT = *titer;
		if( pT->m_key == "du" )
		{
			std::string line= strutil::trim( pT->m_value, "()");
			m_du = strutil::parseString<double>( line );
			break;
		}
	}
};


//write holomorphic 1-form trait m_duv to the string "duv"

inline void CHoloQuadFormEdge::_to_string()
{
	CParser parser( m_string );
	parser._removeToken( "duv" );

	parser._toString( m_string );
	
	std::string line;
	std::stringstream iss(line);
	iss << "duv=(" << m_duv[0] << " " << m_duv[1] << ")";

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

/*! \brief CHoloFormHalfEdge class
 *
 *  HalfEdge class for computing holomorphic 1-forms
 *  Add corner angle m_angle trait
 */

class CHoloQuadFormHalfEdge : public  CHalfEdge
  {
  public:
    /*! corner angle */
	double &   angle() { return m_angle; };

  public:
	/*! CHoloFormHalfEdge constructor */
	 CHoloQuadFormHalfEdge() { m_angle = 0; };
	 /*! CHoloFormHalfEdge destructor */
    ~CHoloQuadFormHalfEdge(){};
	
  protected:
	  /*! corner angle */
	  double   m_angle;

};

/*! \brief CHolomorphicFormMesh class
 *
 *  Mesh class for computing holomorphic forms
 */
template<typename V, typename E, typename F, typename H>
class CHolomorphicQuadraticFormMesh : public CBaseMesh<V,E,F,H>
{
public:
	typedef V CVertex;
	typedef E CEdge;
	typedef F CFace;
	typedef H CHalfEdge;

	typedef CBoundary<V,E,F,H>            CBoundary;
	typedef CLoop<V,E,F,H>                CLoop;
	typedef MeshVertexIterator<V,E,F,H>   MeshVertexIterator;
	typedef MeshFaceIterator<V,E,F,H>	  MeshFaceIterator;
	typedef MeshEdgeIterator<V,E,F,H>     MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef FaceHalfedgeIterator<V,E,F,H> FaceHalfedgeIterator;
	typedef FaceVertexIterator<V,E,F,H>   FaceVertexIterator;
	typedef FaceEdgeIterator<V,E,F,H>	  FaceEdgeIterator;
};


typedef CHolomorphicQuadraticFormMesh<CHoloQuadFormVertex, CHoloQuadFormEdge, CHoloQuadFormFace, CHoloQuadFormHalfEdge> CHoloQuadFormMesh;

} //namespace Holomorphy

unsigned long long Holomorphy::CHoloQuadFormMesh::m_input_traits  = 0;
unsigned long long Holomorphy::CHoloQuadFormMesh::m_output_traits = 0;

} //namespace MeshLib

#endif  _HOLOMORPHIC_QUADRATIC_FORM_MESH_H_