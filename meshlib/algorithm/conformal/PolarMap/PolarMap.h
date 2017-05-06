/*!
*      \file PolarMap.h
*      \brief Algorithm for computing exponential map \f$ z\to e^{z}\f$
*	   \author David Gu
*      \date Documented 10/12/2010
*
*		Computing the exponential map 
*       \f[
*			f(z) = e^{z}
*		\f]
*
*/

/********************************************************************************************************************************
*
*      Polar Map  Class
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       Integration Class
* 
*       David Gu June 27, 2008,  gu@cs.stonybrook.edu
*
*
*      Input:
*         
*           Close Mesh, Open Mesh with UV
*
*      Output:
*
*           Exponential Mapping of the mesh
*
*********************************************************************************************************************************/

/*---------------------------------------------------------------------------------------------------------------------------------
#include "PolarMap/PolarMap.h"

using namespace MeshLib;



int main( int argc, char * argv[] )
{
	CMesh mesh;
	mesh.read_m( argv[1] );

	CMesh open_mesh;
	open_mesh.read_m( argv[2] );

	CPolarMap map( &mesh, &open_mesh );
	map._exponential_map();

	mesh.write_m( argv[3] );

	return 0;
}
----------------------------------------------------------------------------------------------------------------------------------*/

#ifndef _POLAR_MAP_H_
#define _POLAR_MAP_H_


#include  <math.h>
#include <queue>
#include <list>
#include <vector>
#include "PolarMapMesh.h"

namespace MeshLib
{

namespace Holomorphy
{

/*! \brief CPolarMap class
*
*	Algorithm for computing complex exponential map \f$ z\to e^{z} \f$.
*/
template<typename M>
class CPolarMap
{
public:
	/*! CPolarMap constructor
	*
	*  \param pClosedMesh the output closed mesh
	*  \param pOpenMesh   the input open mesh, with vertex uv trait
	*/
	CPolarMap( M * pClosedMesh, M * pOpenMesh );
	/*! CPolarMap destructor */
	~CPolarMap();
	/*! Compute the exponential map*/
	void _exponential_map();

private:
	
	/*! the output closed mesh */
	M * m_pClosedMesh;
	/*! the input open mesh, with vertex uv trait */
	M * m_pOpenMesh;

};


//CPolarMap constructor
//\param pClosedMesh the output mesh
//\param pOpenMesh   the input mesh with vertex uv
template<typename M>
CPolarMap<M>::CPolarMap( M * pClosedMesh, M * pOpenMesh )
{
	m_pClosedMesh = pClosedMesh;
	m_pOpenMesh   = pOpenMesh;
}
//CPolarMap destructor
template<typename M>
CPolarMap<M>::~CPolarMap()
{
}

//Computing the exponential map, 
//Assumption, the range of y value of input vertex uv is from 0 to 1.
template<typename M>
void CPolarMap<M>::_exponential_map()
{
	/*
		for(CPMMesh::MeshVertexIterator viter( m_pOpenMesh ); !viter.end(); viter ++ )
		{
			CPolarMapVertex * pV = *viter;
			CPoint2 p = pV->uv();
		}
	*/

	double x_max = -1e+10;
	for( M::MeshVertexIterator viter( m_pOpenMesh ); !viter.end(); viter ++ )
	 {
		 M::CVertex * pV = *viter;
		 CPoint2 p = pV->uv();
		 x_max = (x_max > p[0])?x_max:p[0];
	}

	printf("X max is %f\n", x_max );

	for( M::MeshVertexIterator viter( m_pOpenMesh ); !viter.end(); viter ++ )
	 {
		 M::CVertex * pV = *viter;
		 CPoint2 p = pV->uv();
		 p[0] -= x_max;
		 pV->uv() = p;
	}

	for( M::MeshVertexIterator viter( m_pOpenMesh ); !viter.end(); viter ++ )
	 {
		 M::CVertex * pV = *viter;
		 CPoint2 p = pV->uv();


		double ang = p[1] * 2 * PI;
		double rad = exp( p[0]* 2 * PI );

		 CPoint2 uv( rad * cos(ang), rad * sin(ang ));
		 pV->uv() = uv;

	}

	

	for( M::MeshVertexIterator viter( m_pOpenMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		int id = pV->father();
		M::CVertex * pW = m_pClosedMesh->idVertex( id );
		pW->uv() = pV->uv();
	}
	

};

} //namespace Holomorphy
} //namespace MeshLib
#endif