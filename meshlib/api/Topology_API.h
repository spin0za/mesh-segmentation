#ifndef _TOPOLOGY_API_H_

/*!
 *	Trait IO 
 */
#include "parser/traits_io.h"

/*************************************************************************************************************************
*
*	Topological methods
*
*************************************************************************************************************************/

/*! Integrating holomorphic 1-form on a simply connected surface
 */
#include "Topology/Integration/Integration.h" //Integrating holomorphic 1-form on a simply connected surface
/*! Integrating combinatorial 1-form on a simply connected surface
 *  to obtain a square tiling
 */
#include "Topology/Integration/SquareTilingIntegrationMesh.h" //Integrating combinatorial 1-form Mesh
#include "Topology/Integration/SquareTilingIntegration.h"     //Integrating combinatorial 1-form on a simply connected surface to obtain a square tiling
/*! Cut Graph
*/
#include "Topology/CutGraph/CutGraph.h" //Cut Graph
/*!	Slice the open along sharp edges to form another mesh - Wedge mesh
 */
#include "Topology/Wedge/WMesh.h" //Slice the open along sharp edges to form another mesh - Wedge mesh

/*!	Homology Basis
 */
#include "Topology/CutGraph/Homology.h" //Homology
/*!	Cohomology Group basis , closed one form
 */
#include "Topology/Cohomology/Cohomology.h" 
/*!	Cohomology Group basis for multiply connected domain, closed one form
 */
#include "Topology/Cohomology/DomainCohomology.h" 
/*! Puncture a hole in the center of the mesh
 */
#include "Topology/puncture/puncture.h" //Puncture a hole in the center of the mesh
/*!
 *	Double covering, double cover a mesh with boundaries to a closed symmetric mesh
 */
#include "Topology/DoubleCovering/DoubleCoveringMesh.h"

/*!	Shortest Path
 */
#include "Riemannian/ShortestPath/ShortestPath.h" //shortest path

/*! Linear combination of 1-forms
*/
#include "Topology/Integration/LinearCombination.h"
/*! Linear combination of holomorphic 1-forms
*/
#include "Topology/Integration/LinearCombinationOneForm.h"

/*! Stitch two meshes along their common boundaries, Fill the center hole of a planar mesh
 */
#include "Topology/Stitch/Stitch.h" //Stitch two meshes along their common boundaries

/*! Dijkstra, Compute the tight fundamental group generators
*/
#include "Riemannian/Dijkstra/Dijkstra.h"
/*!
 *	harmonic form mesh
 */
#include "Conformal/ClosedForm/HarmonicClosedFormMesh.h"


/************************************************************************************************************************************
*
*	Topological Algorithms
*
************************************************************************************************************************************/

/*!	Compute the cut graph of a mesh
 *
 */
void _cut_graph( const char * _input, const char * _output );

/*!	compute the shortest path connecting an inner boundary to the exterior boundary
 *
 */
void _cut_domain( const char * _domain_mesh, const char * _mesh_with_cut );

/*! Slice the mesh open along the sharp edges
 *
 */
void _slice( const char * _closed_mesh, const char * _open_mesh );

/*!	compute the homology group basis
 *
 */
void _homology( const char * _input, const char * _output );

/*!	compute the closed one form basis
 *
 */
void _cohomology_one_form( const char * _closed_mesh, const char * _open_mesh, const char * _closed_1_form, const char * _integrated );

/*!	compute the closed one form basis for multiply connected domain
 *
 */
void _cohomology_one_form_domain( const char * _input_mesh, const char * _sliced_mesh, const char * _closed_1_form, const char * _integrated_mesh );

/*! Puncture a hole in the center of the mesh
 *
 */
void _puncture( const char * _input_mesh, const char * _mesh_with_hole );

/*! Fill the hole in the center
 *
 */
void _fill_puncture( const char * _mesh_with_hole, const char * _filled_mesh );


/*! Koebe's Iteration Method: Fill the central hole
 */
void Koebe_fill_centeral_hole( const char * _mesh_with_hole, const char * _filled_mesh );

/*! Koebe's iteration Method: Remove faces with a given segment ID
 */
void Koebe_remove_segment( const char * _input_mesh, const char * segment_id, const char * _output_mesh );

/*! Koebe's iteration method: Copy UV, from UV_Mesh to Position_Mesh
 */
void Koebe_copy_uv( const char * _uv_mesh, const char * _pos_mesh, const char * _output_mesh );

/*! Compute the double covering mesh
 */
void _double_covering( const char * _open_mesh, const char * _doubled_mesh );

/*! Compute the shortest basis using Erickson's method
 */
void _shortest_homology_basis( const char * _input, const char * _vertex_id, const char * _basis_mesh );

/*!	Linear combination of holomorphic 1-forms
*
*/
void _linear_combination( int argc, char * argv[] );
/*! Linear combination of 1-forms
 *
 */
void _linear_combination_one_form( int argc, char * argv[] );

/*!	Integrate a holomorphic 1-form on a simply connected surface
 *
 */
void _integration( const char * _form_input, const char * _domain_input, const char *_output );
/*!	Integrate a combinatorial harmonic 1-form on a simply connected surface
 *  to obtain a square tiling 
 */
void _square_tiling_integration( const char * _form_input, const char * _domain_input, const char *_output );

#endif _TOPOLOGY_API_H_
