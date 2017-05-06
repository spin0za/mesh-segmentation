/*! \file HolomorphicFormMesh.h
*    \brief Mesh for Computing Holomorphic Forms
*	 \author David Gu
*    \date Document 10/12/2010
*/

#ifndef  _HOLOMORPHIC_FORM_MESH_H_
#define  _HOLOMORPHIC_FORM_MESH_H_

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

	Holomorphic Form Face Trait


-------------------------------------------------------------------------------------------------------------------------------------------*/

/*!
 *	\brief CHoloFormVertex class
 *
 *	Vertex class for computing holomorphic forms
 *  add vertex huv and mu
 */

class CHoloFormVertex: public CVertex
{
public:
	/*! CHoloFormVertex constructor */
	CHoloFormVertex() {};
	/*! CHoloFormVertex destructor  */
	~CHoloFormVertex() {};
	
	/*! conformal parameter */
	std::complex<double> & z() { return m_z; };
	/*! vertex Beltrami coefficient */
	std::complex<double> &mu() { return m_mu; };

protected:
	/*! conformal parameter */
	std::complex<double>  m_z;
	/*! vertex mu */
	std::complex<double> m_mu;

};


/*!
 *	\brief CHoloFormFace class
 *
 *	Face class for computing holomorphic forms
 *  add face area, vector representation of closed 1-form and face normal traits
 */

class CHoloFormFace: public CFace
{
protected:
	/*! face area trait */
	double  m_area;
	/*! vector representation of closed 1-form */
	CPoint  m_du;
	/*! normal to the face */
	CPoint  m_normal;

public:
	/*! CHoloFormFace constructor */
	CHoloFormFace() { m_area = 0; };
	/*! CHoloFormFace destructor  */
	~CHoloFormFace() {};

public:
	/*! Face area */
	double & area()   { return m_area;   };
	/*! normal vector to the face */
	CPoint & normal() { return m_normal; };
	/*! Vector representation of a closed 1-form on the face */
	CPoint & du()     { return m_du;     };
};


/*-----------------------------------------------------------------------------------------------------------------------------------------

	Holomorphic Form Edge Trait


-------------------------------------------------------------------------------------------------------------------------------------------*/

 /*! \brief CHoloFormEdge class
  *
  *  Edge class for computing holomorphic forms
  *  add harmonic 1-form trait m_du
  *  holomorphic 1-form trait m_duv
  *  edge length m_length traits.
  */
class CHoloFormEdge : public  CEdge
{
  protected:
	  /*! edge length */
	  double  m_length;
	  /*! harmonic 1-form */
	  double  m_du;
	  /*! holomorphic 1-form */
	  CPoint2 m_duv;
  
  public:
	  /*! CHoloFormEdge constructor */
      CHoloFormEdge() { m_length = 0.0; m_du = 0;  };
	  /*! CHoloFormEdge destructor */
      ~CHoloFormEdge(){};
	  /*! length of the edge */
	  double  & length() { return m_length; };
	  /*! harmonic 1-form */
	  double  & du()     { return m_du;     };
	  /*! holomorphic 1-form */
	  CPoint2 & duv()    { return m_duv;    };
	  /*! read harmonic 1-form from edge string */
	  void _from_string();
	  /*! write holomorphic 1-form to edge string */
	  void _to_string();
};

//read harmonic 1-form trait "du" to the trait m_du
inline void CHoloFormEdge::_from_string()
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

inline void CHoloFormEdge::_to_string()
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

class CHoloFormHalfEdge : public  CHalfEdge
  {
  public:
    /*! corner angle */
	double &   angle() { return m_angle; };

  public:
	/*! CHoloFormHalfEdge constructor */
	 CHoloFormHalfEdge() { m_angle = 0; };
	 /*! CHoloFormHalfEdge destructor */
    ~CHoloFormHalfEdge(){};
	
  protected:
	  /*! corner angle */
	  double   m_angle;

};

/*! \brief CHolomorphicFormMesh class
 *
 *  Mesh class for computing holomorphic forms
 */
template<typename V, typename E, typename F, typename H>
class CHolomorphicFormMesh : public CBaseMesh<V,E,F,H>
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
};


typedef CHolomorphicFormMesh<CHoloFormVertex, CHoloFormEdge, CHoloFormFace, CHoloFormHalfEdge> CHoloFormMesh;

} //namespace Holomorphy

unsigned long long Holomorphy::CHoloFormMesh::m_input_traits  = EDGE_DU;
unsigned long long Holomorphy::CHoloFormMesh::m_output_traits = EDGE_DUV;

} //namespace MeshLib
#endif  _HOLOMORPHIC_FORM_MESH_H_