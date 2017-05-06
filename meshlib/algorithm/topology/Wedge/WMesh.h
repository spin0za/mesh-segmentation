/*!
*      \file WMesh.h
*      \brief Algorithm for slicing a mesh along sharp edges to form a new mesh
*	   \author David Gu
*      \date Document 10/11/2010
*
*		Slice a mesh along sharp edges to form a new mesh. Input is a mesh with sharp edges,
*		the output is the open mesh sliced along the sharp edges. 
*/


/*******************************************************************************
*      Wedge Mesh  Class
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       Wedge Mesh Class
* 
*       David Gu June 27, 2008, gu@cs.stonybrook.edu
*   
*    Input:
*       
*        A mesh with sharp edges 
*
*    Output
*		
*        The mesh is sliced along the sharp edges to form an open mesh
*
*******************************************************************************/


/*------------------------------------------------------------------------------------------------------------------------------------

//Example Code

#include "mesh/mesh.h"
#include "Wedge/WMesh.h"

using namespace MeshLib;

int main( int argc, char * argv[] )
{
	CMesh mesh;
	mesh.read_m( argv[1] );

	CWMesh wmesh( & mesh );
    wmesh.Slice();

    wmesh.wmesh()->write_m( argv[2] );

  return 0;
}

------------------------------------------------------------------------------------------------------------------------------------*/

#ifndef  _WMesh_H_
#define  _WMesh_H_


#include <assert.h>
#include <math.h>

#include <iostream>
#include <map>
#include <list>

#include "WedgeMesh.h"


namespace MeshLib{

namespace Topology
{
/*!
 *  \brief CWedge class
 *
 *   A wedge is a union of corners sharing the same apex vertex. The corners attaching to the same vertex are
 *   partitioned to different wedges by sharp edges or boundary edges.
 */
class CWedge  
{
public:
	/*! CWedge constructure */
	CWedge(){};
	/*! CWedge destructor */
	~CWedge(){};

	/*! CWedge constructor
	*	\param pMesh the input mesh
	*	\param half_edge the first corner in this wedge
	*/
	CWedge( CSMesh * pMesh, CWedgeHalfEdge *half_edge);//constructor 
	/*!	The vertex of the current wedge
	*/
	CWedgeVertex *  vertex() { return m_vertex; };						//the vertex of the wedge
	/*!	The most clockwise corner of the current wedge
	 */
	CWedgeHalfEdge * mostClwHalfEdge(){ return m_mostClwHEdge;};			//return the most Clw HEdge of this wedge
	/*!	The most counter clockwise corner of the current wedge
	 */
	CWedgeHalfEdge * mostCcwHalfEdge(){return m_mostCcwHEdge;};			//return the most Ccw HEdge of this wedge
	/*!	Each wedge corresponds to a vertex in wmesh, the associated wmesh vertex of current wedge
	 */
	CWedgeVertex * & wvertex() { return m_wvertex; };					//the new vertex on the wedge solid

protected:
	/*!	The current wedge mesh
	*/
	CSMesh				  * m_mesh;
	/*!	vertex of all the corners in the wedge
	*/
	CWedgeVertex	      * m_vertex;
	/*!	The most CCW corner of the wedge.
	 */
	CWedgeHalfEdge		  * m_mostCcwHEdge;
	/*! The most CLW corner of the wedge.
	 */
	CWedgeHalfEdge        * m_mostClwHEdge;
	/*!	New vertex on the wedge mesh, corresponding to this wedge.
	 */
	CWedgeVertex	      * m_wvertex;	//new vertex on wedge solid

};

//Wedge Solid, Support mesh slicing
/*!
 *	\brief CWMesh class
 *
 *  Converting a mesh with sharp edges to a new mesh, such that 
 *  1. each wedge becomes a new vertex
 *  2. each face of the old mesh becomes a new face in the following way: each corner of the old face belongs to an wedge,
 *     the three wedges are connected to a new face.
 */

class CWMesh : public CSMesh
{

public:
	/*! CWMesh constructor */
	//constructor and destructor
	CWMesh();
	/*! CWMesh constructor
	* \param pMesh the input mesh with sharp edges
	*/
	CWMesh( CSMesh * pMesh );
	/*! CWMesh destructor */
	~CWMesh();
	/*! The newly constructed mesh */
	CSMesh * wmesh() { return &m_wmesh; };	
	/*! Slice the input mesh along the sharp edges. */
	void   Slice();

private:
	/*! Construct a halfedge structure, each vertex is a wedge. */
	void		_construct();
	/*! Convert the halfedge structure with wedge vertices to a common mesh. */
	void		_convert();
	/*! The topological valence of the vertex, number of sharp edges or boundary edges,
	 *  \param vertex input vertex
	 */
	int		    __topovalence( CWedgeVertex * vertex );			//compute topological valence of vertex

	/*! The input mesh. */
	CSMesh *		m_pMesh;
	/*! The output converted mesh. */
	CSMesh		    m_wmesh;
	/*! list of wedges. */
	std::list<CWedge*> m_wedges;							//buffer for wedges
};



}//name space Topology
}//name space MeshLib

#endif //_MESHLIB_SOLID_H_ defined



