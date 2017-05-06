/*!
*      \file ShortestPath.h
*      \brief Algorithm for shortest paths
*	   \author David Gu
*      \date Document 10/07/2010
*
*		Compute the shortest path from one inner boundary loop to the exterior boundary loop. The input is a mesh with 
*       multiple boundary components; the output are meshes with shortest paths, labeled as sharp edges,  between boundary 
*       k to boundary 0 (exterior boundary).
*/

/*******************************************************************************
*      Shortest Path between two boundary loops
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       Compute the shortest path from one inner boundary loop to the exterior boundary loop
* 
*       David Gu June 27, 2008, gu@cs.stonybrook.edu
*
*	Input
*       A mesh with multiple boundaries
*	Output
*       The shortest path labeld as sharp edges, between boundary k to boundary 0 (exterior boundary).
*
*******************************************************************************/

/*-------------------------------------------------------------------------------------------------------------------------------

#include <math.h>
#include "mesh/mesh.h"
#include "ShortestPath/ShortestPath.h"

using namespace MeshLib;

int main( int argc, char * argv[] )
{
 
  if( strcmp( argv[1], "-cut" ) == 0 )
  {

	CMesh cmesh;
	cmesh.read_m( argv[2] );

	CShortestPath cut( & cmesh );
	cut._cut( argv[3] );

  return 0;
  }


}

--------------------------------------------------------------------------------------------------------------------------------*/
#ifndef _SHORTEST_PATH_H_
#define _SHORTEST_PATH_H_

#include  <math.h>
#include <queue>
#include "Mesh/boundary.h"
#include "Mesh/iterators.h"
#include "ShortestPathMesh.h"

namespace MeshLib
{
/*!
 * \brief CShortestPath class
 * 
 * Compute the shortest path between one interior boundary component and the exterior boundary component
 * 
 */
  class CShortestPath
  {
  public:
    /*!
	 * CShortestPath default constructor
	 * \param pMesh the input mesh
	 */
    CShortestPath( CSPMesh * pMesh );
	/*! CShortestPath destructor
	 *
	 */
    ~CShortestPath();
	/*!	Compute the shortest path from boundary loop k to boundary loop 0 (exterior boundary loop ),
	 *  the path is labeled as sharp edges and output to "output_name_k.m"
	 *  \param prefix the prefix of output mesh name
	 */
    void _cut( const char * prefix );
  
  protected:
    /*! Pointer to the input mesh
	 */
    CSPMesh * m_pMesh;
	/*!	Boundary loops of the input mesh
	 */
	CSPMesh::CBoundary m_boundary;
	/*! Edges on the shortest path 
	 */
	std::list<CSPEdge*> m_cuts;
	/*!	Compute the shortest path between source loop to the target loop
	 * \param source the starting boundary loop
	 * \param target the ending boundary loop
	 */
	void _trace(CSPMesh::CLoop * source, CSPMesh::CLoop * target);
  };
}
#endif