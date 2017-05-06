/*!
*      \file Stitch.h
*      \brief Algorithm for Stiching two meshes
*	   \author David Gu
*      \date Documented on 10/12/2010
*
*	   Algorithm for stitching two meshes along their common boundaries.
*/

/*******************************************************************************
*      Stich mesh patches together
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       Locate a boundary loop, generate a mesh patch, stich the mesh patch with the original mesh
* 
*       David Gu June 27, 2008, gu@cs.stonybrook.edu
*
*	Input
*       A mesh with multiple boundaries
*	Output
*       Stich mesh patches with the mesh to form the new mesh
*
*******************************************************************************/

/*-------------------------------------------------------------------------------------------------------------------------------

#include <math.h>
#include "Stitch/Stitch.h"

using namespace MeshLib;

int main( int argc, char * argv[] )
{
 
	//for stitching purpose
	if( strcmp( argv[1], "-hole2poly" ) == 0 )
	{
		CSTMesh mesh;
		mesh.read_m( argv[2] );
		CStitch stitch( &mesh );
		stitch._hole_2_poly( argv[3] );
		return 0;
	}

	if( strcmp( argv[1], "-off2mesh" ) == 0 )
	{
		CTopoMesh mesh;
		mesh.read_off( argv[2] );
		mesh.write_m( argv[3] );
		return 0;
	}

	if( strcmp( argv[1], "-merge" ) == 0 )
	{

		CSTMesh mesh;
		mesh.read_m( argv[2] );

		CSTMesh patch;
		patch.read_m( argv[3] );
	
		CStitch stitch( &mesh );
		stitch._stitch( patch, argv[4] );
		return 0;
	}
	return 0;

}

--------------------------------------------------------------------------------------------------------------------------------*/
#ifndef _STITCH_H_
#define _STITCH_H_

#include "StitchMesh.h"
#include "Riemannian/Delaunay/DelaunayMesh.h"

namespace MeshLib
{
namespace Topology
{
	/*! \brief CStitch class
	*
	*	Algorithm for stitching two meshes along their two common boundaries.
	*/
  class CStitch
  {
  public:
	  /*! CStich constructor
	  *  \param pMesh input base mesh
	  */
    CStitch( CSTMesh * pMesh );
	/*! CStitch destructor */
    ~CStitch();

	/*! Converting a boundary loop to the poly file
	* \param output_name output poly file name 
	*/
	void _hole_2_poly( const char * output_name );
	/*! Stitch mesh patch to the current base mesh 
	 * \param patch the input patch mesh
	 * \param output_name the output mesh file name
	 */
	void _stitch( CSTMesh & patch, const char * output_name );
	
	/*! Fill the biggest hole
	 */
	void _fill_the_biggest_hole(  const char * output_name );

	/*! Punch a hole by removing all faces with segment id 
	 *  \param segment_id, the segment intended to be removed
	 *  \param output_mesh name
	 */
	void _remove_segment(  int segment_id, const char * output_name );

	/*! Copy uv
	 *  \param CSTMesh mesh
	 */
	void _copy_uv( CSTMesh * mesh );

  protected:
	/*! the base mesh */
	CSTMesh * m_pMesh;
	/*! the boundary of the base mesh */
	CSTMesh::CBoundary m_boundary;
	/*! union the current mesh with patch, the output is union_mesh 
	 *	
	 *	\param patch the patch mesh
	 *  \param union_mesh the union mesh
	 */
	void _union( CDTMesh * patch, CSTMesh & union_mesh );

  private:

	  //planar area circled by a loop
	  /*! planar area circled by a loop
	  *   \param hl the boundary loop
	  *   \return the area of the domain circled by hl
	  */
	  double _area( std::list<CStitchHalfEdge*>  & hl );

  };
} //namespace Topology
} //namespace MeshLib
#endif