/*!
*      \file SphericalHarmonicMap.h
*      \brief Algorithm for Spherical Harmonic Maps
*	   \author David Gu
*      \date Document 12/15/2010
*
*		Computing spherical harmonic maps
*
*/

/********************************************************************************************************************************
*
*      Spherical Harmonic Map Class
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       Compute Spherical Harmonic Maps
* 
*       David Gu December 15, 2010,  gu@cs.stonybrook.edu
*
*
*      Input:
*         
*           Genus zero closed surface
*
*      Output:
*
*           Image of the spherical harmonic mapping
*
*********************************************************************************************************************************/

/*---------------------------------------------------------------------------------------------------------------------------------
#include <math.h>
#include "mesh/mesh.h"
#include "SphericalHarmonicMap/SphericalHarmonicMap.h"

using namespace MeshLib;

int main( int argc, char * argv[] )
{

  if( strcmp( argv[1], "-spherical_harmonic_map" ) == 0 )
  {
	CMesh mesh;
	mesh.read_m( argv[2] ); //original mesh

	CSphericalHarmonicMap map( &mesh );
	map.map();

	mesh.write_m( argv[3] );
	return 0;
  }

}
----------------------------------------------------------------------------------------------------------------------------------*/

#ifndef _SPHERICAL_HARMONIC_FORM_H_
#define _SPHERICAL_HARMONIC_FORM_H_

#include  <math.h>
#include <queue>
#include "Structure/Structure.h"
#include "SphericalHarmonicMapMesh.h"

namespace MeshLib
{

namespace Holomorphy
{

	/*! \brief CSphericalHarmonicMap class
	*  
	*  Algorithm for computing spherical harmonic maps
	*/
  template<typename M>
  class CSphericalHarmonicMap
  {
  public:
	  /*! CSphericalHarmonicMap constructor
	  *
	  *	\param pMesh input genus zero closed mesh
	  */

    CSphericalHarmonicMap( M * pMesh  );
	/*!
	 * CSphericalHarmonicMap destructor
	 */
    ~CSphericalHarmonicMap();
	/*!
	 *	Compute the spherical Harmonic Map, calling function
	 */
	void map();
	/*!
	 *	compute spherical harmonic map
	 *  \param step_length step length
	 *  \param threshold   threshold
	 */
    void _harmonic_map(double step_length, double threshold );
	/*!
	 *	compute Tutte map
	 *  \param step_length step length
	 *  \param threshold   threshold
	 */

	void _Tutte_harmonic_map( double step_length, double threshold );

  protected:
    /*!
	 *	Input genus zero closed mesh
	 */ 
    M * m_pMesh;


	/*!
	 *	Compute the Tutte edge weight
	 */
	void _edge_Tutte_weight();
	/*!
	 *	Compute the edge weight
	 */
	void _edge_weight();

	/*!
	 *	Calculate the normals at each vertex
	 */
	void _calculate_normal();

	/*!
	 *	Calculate the Gauss map
	 */
	void _calculate_Gauss_map();

	/*!
	 *	Calculate the harmonic energy
	 */
	double _harmonic_energy();
	/*!
	 *	Calculate the Laplacian at each vertex
	 */
	void   _calculate_Laplacian();

	/*!
	 *	Normalization
	 */
	void _normalize();

	/*!
	*	Projecting the Laplacian to the tangent space
	*/
	void _projection();
	/*!
	 *	Update the mapping
	 *  \param step_length step length
	 */
	void _update( double step_length );


  };

template<typename M>
CSphericalHarmonicMap<M>::CSphericalHarmonicMap( M * pMesh )
{
	m_pMesh = pMesh;

};

//Destructor
template<typename M>
CSphericalHarmonicMap<M>::~CSphericalHarmonicMap()
{
};

//calculate the edge Tutte weight
template<typename M>
void CSphericalHarmonicMap<M>::_edge_Tutte_weight()
{
  for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end();  eiter ++ )
  {
	  M::CEdge * e = *eiter;
	  e->weight() = 1.0;
  }
}

//calculate the edge weight
template<typename M>
void CSphericalHarmonicMap<M>::_edge_weight()
{
	COperator<M> pC( m_pMesh );
	pC._embedding_2_metric();
	pC._metric_2_angle();
	pC._angle_2_Laplace();	//compute edge weight
}

//Calculate Vertex Normal
template<typename M>
void CSphericalHarmonicMap<M>::_calculate_normal()
{
	for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); ++ fiter )
	{
		M::CFace * f = *fiter;

		M::CHalfEdge * he = m_pMesh->faceHalfedge( f );
		CPoint p[3];
		
		for( int i = 0; i  < 3; i ++ )
		{
			M::CVertex * v = m_pMesh->halfedgeTarget( he );
			p[i] = v->point();
			he = m_pMesh->faceNextCcwHalfEdge( he );
		}
		
		CPoint n = (p[1]-p[0])^(p[2]-p[0]);

		f->area() = n.norm();
		n/=n.norm();
		f->normal() = n;
	}

	for( M::MeshVertexIterator viter(m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * v = *viter;
		CPoint n(0,0,0);
		double w = 0;
		for( M::VertexFaceIterator vfiter( v ); !vfiter.end() ; ++ vfiter )
		{
			M::CFace * f = *vfiter;
			n = n + f->normal() * f->area();
			w = w + f->area();
		}
		n/=w;
		n/=n.norm();
		v->normal() = n;
	}
}

//Calculate the Gauss map of the surface
template<typename M>
void CSphericalHarmonicMap<M>::_calculate_Gauss_map()
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * v = *viter;
		v->u() = v->normal( );
		assert( fabs( v->u().norm() - 1 ) < 1e-4 );
	}
}

//Calculate Harmonic Energy
template<typename M>
double CSphericalHarmonicMap<M>::_harmonic_energy()
{
	double sum = 0;

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); ++ eiter )
	{
		M::CEdge * e  = *eiter;
		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );
		CPoint d = v1->u() - v2->u();
		sum = sum + e->weight() * ( d * d );
	}

	return sum;
}

//calculate the Laplacian
template<typename M>
void CSphericalHarmonicMap<M>::_calculate_Laplacian()
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * v = *viter;
		CPoint    L(0,0,0);
		double sw = 0;

		for( M::VertexOutHalfedgeIterator vhiter( m_pMesh, v ); !vhiter.end(); ++ vhiter )
		{
			M::CHalfEdge * he = *vhiter;
			M::CVertex * w = m_pMesh->halfedgeTarget( he );
			assert( m_pMesh->halfedgeSource(he) == v );
			M::CEdge * e = m_pMesh->halfedgeEdge( he );
			double wg = e->weight();
			L = L + ( w->u() - v->u() ) * wg;
			sw += wg;
		}
		L/= sw;
		v->L() = L;
	}
}

//normalization of the map
template<typename M>
void CSphericalHarmonicMap<M>::_normalize()
{
	int count = 0;

	CPoint center(0,0,0);
	
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * v = *viter;
		center += v->u();
		count ++;
	}

	center /= (double) count;

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * v = *viter;
		CPoint u = v->u();
		u = u - center;
		u/=u.norm();
		v->u() = u;
	}
}

//update the mapping
template<typename M>
void CSphericalHarmonicMap<M>::_update( double step_length )
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * v = *viter;
		CPoint L = v->L();
		L = L * step_length;
		v->u( ) += L;
		v->u( ) /= v->u().norm();
	}
}

//project the Laplacian to the tangential space
template<typename M>
void CSphericalHarmonicMap<M>::_projection( )
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * v = *viter;
		CPoint L = v->L();
		CPoint u = v->u();
		u/= u.norm();
		double dot = L * u;
		L= L -  u * dot;
		dot = L * u;
		v->L() = L;
	}
}

//Tuette Harmonic Mapping
template<typename M>
void CSphericalHarmonicMap<M>::_Tutte_harmonic_map( double step_length, double threshold )
{
	_calculate_normal();
	_calculate_Gauss_map();

	_edge_Tutte_weight();

	double e0 = _harmonic_energy();
	int count = 0;
	while( true )
	{
		_calculate_Laplacian();
		_projection();
		_update( step_length );
		_normalize();

		double e1 = _harmonic_energy();
		if( count ++ % 100 == 0 )
			printf("%f\n", e1 );

		if( fabs( e0 - e1 ) < threshold ) 
			break;
		e0 = e1;
	}
//	for( MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
//	{
//		CVertex * v = *viter;
//		v->point() = v_u( v );
//	}

}

//compute the spherical harmonic map
template<typename M>
void CSphericalHarmonicMap<M>::_harmonic_map( double step_length, double threshold )
{
/*	_calculate_norm();
	_calculate_Gauss_map();
*/
	//compute cotangent edge weight
	_edge_weight();

	double e0 = _harmonic_energy();

	while( true )
	{
		_calculate_Laplacian();
		_projection();
		_update( step_length );
		_normalize();

		double e1 = _harmonic_energy();
		printf("%f\n", e1 );

		if( fabs( e0 - e1 ) < threshold ) 
			break;
		e0 = e1;
	}
}

//compute the spherical harmonic map, calling function
template<typename M>
void CSphericalHarmonicMap<M>::map()
{
	printf("Tuette Map\n");
	_Tutte_harmonic_map( 9e-1, 0.001 );
	printf("Harmonic Map\n");
	_harmonic_map( 5e-1,1e-3);
}




} //namespace Holomorphy

} //namespace MeshLib

#endif _SPHERICAL_HARMONIC_MAP_H_