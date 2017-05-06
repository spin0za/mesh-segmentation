#ifndef _DELAUNAY_API_H_
#define _DELAUNAY_API_H_

/*!
 *	Trait IO 
 */
#include "parser/traits_io.h"

/************************************************************************************************************************************
*
*	Delaunay Meshing
*
************************************************************************************************************************************/

#include "Riemannian/Delaunay/DelaunayMesh.h"
#include "Riemannian/Delaunay/DelaunayRemesh.h"
#include "Riemannian/Delaunay/DelaunayRefineMesh.h"
#include "Riemannian/Delaunay/QuadTreeDelaunayRemesh.h"

/*! Planar Mesh Generation
*
*/

void _Rupert_Delaunay_Refinement( const char * poly_file, const char * mesh_file );
/*! Surface Remeshing based on Delaunay Refinement Method
*
*/
void _Surface_Delaunay_Remesh( const char * _input_mesh, const char * _output_mesh );
/*! Surface meshing Delaunay refinement, only insert new vertices, no vertex is removed
*
*/
void _Surface_Delaunay_Mesh_Refinement( const char * _input_file, const char * _output_file );
/*!
 *	Surface Delaunay remshing using Quad Tree data structure to improve the robustness
 */
void _Surface_Delaunay_Quad_Tree_Remesh( const char * _input_mesh, const char * _output_mesh );


#endif _DELAUNAY_API_H_
