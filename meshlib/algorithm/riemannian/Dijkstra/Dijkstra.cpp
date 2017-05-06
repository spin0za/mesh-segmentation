/*! \file Dijkstra.cpp
*   \brief Implementation of CDijkstra class
*	\author David Gu
*   \date Documented on 03/20/2011	
*
*	Compute the tight fundamental group generators
*   1. Compute the shortest non-separating, essential loop through a given point
*   2. Cut the surface along the loop
*   3. Compute the shortest non-separating, essential loop, connecting two corresponding points on the boundary.
*/

#include "Dijkstra.h"

using namespace MeshLib;

CDijkstra::CDijkstra( CDKMesh * pMesh, CDKMesh * pDomain ):m_boundary( pDomain )
{
	m_pMesh   = pMesh;
	m_pDomain = pDomain;

};

CDijkstra::~CDijkstra()
{
};





