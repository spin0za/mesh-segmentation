/*!
*      \file Puncture.h
*      \brief Algorithm for punching a hole in the center of the mesh, and the reverse operation: filling the hole.
*	   \author David Gu
*      \date Documented on 10/12/2010
*
*/

/*******************************************************************************
*      Puncture Hole on a mesh/ fill a hole of the mesh
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       Locate the vertex in the center, puncture an adjacent face / locate the smallest boundary, fill it with a face
* 
*       David Gu June 27, 2008, gu@cs.stonybrook.edu
*
*	Input
*       A mesh/ A mesh with a hole with 3 vertices
*	Output
*       The mesh with a face punctured, the face is adjacent to the center vertex / Fill the hole
*
*******************************************************************************/

/*-------------------------------------------------------------------------------------------------------------------------------

#include "puncture/puncture.h"

using namespace MeshLib;

int main( int argc, char * argv[] )
{
  if( strcmp( argv[1], "-puncture" ) == 0 )
  {

	CPMesh cmesh;
	cPmesh.read_m( argv[2] );

	CPuncture cut( & cmesh );
	cut._puncture();
	
	cmesh.write_m( argv[3] );
	return 0;
  }

  if( strcmp( argv[1], "-fill_hole" ) == 0 )
  {

	CPMesh cmesh;
	cmesh.read_m( argv[2] );

	CPuncture cut( & cmesh );
	cut._fill_hole();
	
	cmesh.write_m( argv[3] );
	return 0;
  }

}

--------------------------------------------------------------------------------------------------------------------------------*/
#ifndef _PUNCTURE_H_
#define _PUNCTURE_H_

#include <queue>
#include "PunctureMesh.h"

namespace MeshLib
{

namespace Topology
{
/*! \brief CPuncture class
*
*	Algorithm for punching a hole in the center of the mesh
*   and filling the hole. This is used to compute the Riemann mapping
*   using slitMap.
*/
  template<typename M>
  class CPuncture
  {
  public:
	  /*! CPuncture class constructor
	  * \param pMesh the input mesh
	  */
    CPuncture( M * pMesh );
	/*! CPuncture class destructor */
    ~CPuncture();

	/*! puncture a hole in the center of the mesh */
	void _puncture();
	/*! fill the holed punctured */
	void _fill_hole();
  
  protected:
	  /*! input mesh */
    M * m_pMesh;
	/*! the boundary of the input mesh */
	typename  M::CBoundary m_boundary;

  };

//CPuncture constructor
//\param pMesh the input mesh
template<typename M>
CPuncture<M>::CPuncture( M * pMesh ): m_pMesh( pMesh ), m_boundary( m_pMesh )
{
}

//CPuncture destructor
template<typename M>
CPuncture<M>::~CPuncture()
{
}

//puncture a hole in the center of the mesh
//we use bread first search method, to find the most inner vertex on the mesh
// and treat it as the center of the mesh.
template<typename M>
void CPuncture<M>::_puncture()
{
/*
	if( m_boundary.loops().size() != 1 )
	{
		fprintf(stderr, "Input mesh should be a topological disk\n");
		return;
	}
*/
	std::queue<M::CVertex*> vqueue;

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
    {
		M::CVertex * v = *viter;
	  if( !v->boundary() )
		 v->touched() = false;
	  else
	  {
		  v->touched() = true;
		  vqueue.push( v );
	  }
	
    }
    
	M::CVertex * destiny = NULL;


    while( ! vqueue.empty() )
    {
		M::CVertex * head = vqueue.front();
        vqueue.pop();
		destiny = head;		        
		for( M::VertexVertexIterator vviter( head ); !vviter.end(); ++ vviter )
        {
			M::CVertex * w = *vviter;
            if( w->touched() ) continue;
            w->touched() = true;
            vqueue.push( w );
        }
    }

	printf("Find the center vertex %d\n", destiny->id() );

	M::CFace * punched_face = NULL;

	for( M::VertexFaceIterator vfiter( destiny ); !vfiter.end(); vfiter ++ )
	{
		M::CFace * pF = *vfiter;
		punched_face = pF;
		break;
	}
	assert( punched_face != NULL );
	m_pMesh->deleteFace( punched_face );
}

//filling the small hole
template<typename M>
void CPuncture<M>::_fill_hole()
{
/*
	if( m_boundary.loops().size() != 2 )
	{
		fprintf(stderr, "Input mesh should be a topological disk with a small hole\n");
		return;
	}
*/
	M::CLoop * pL = m_boundary.loops().back();

	std::vector<M::CVertex*> verts;

	for( std::list<M::CHalfEdge*>::iterator hiter = pL->halfedges().begin(); hiter != pL->halfedges().end(); hiter ++ )
	{
		M::CHalfEdge * pH = *hiter;
		verts.push_back( m_pMesh->halfedgeTarget( pH ) );
	}


	int fid = -1;

	for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		M::CFace * pF = *fiter;
		fid = (fid> pF->id())? fid: pF->id();
	}

	for( size_t j = 0; j < verts.size() - 2; j ++ )
	{
		M::CVertex *v[3];
		
		v[0] = verts[0];
		v[2] = verts[j+1];
		v[1] = verts[j+2];
		m_pMesh->createFace( v, ++fid );
	}
}

} //namespace Topology
} //namespace MeshLib
#endif