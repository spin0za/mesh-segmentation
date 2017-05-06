#include "Topology_API.h"

using namespace MeshLib;

/**********************************************************************************************************************************************
*
*	Topological Algorithms
*	
*
**********************************************************************************************************************************************/

/*!	Integrate a holomorphic 1-form on a simply connected surface
 *
 */

void _integration( const char * _form_input, const char * _domain_input, const char *_output )
{
	Topology::CIMesh holo_mesh;
	Topology::CIMesh fund_mesh;

	holo_mesh.read_m( _form_input   );
	fund_mesh.read_m( _domain_input );

	Topology::CIntegration<Topology::CIMesh> integrator( & holo_mesh, & fund_mesh );
	integrator._integrate();

	fund_mesh.write_m( _output );
}

/*!	Integrate a combinatorial harmonic 1-form on a simply connected surface
 *  to obtain a square tiling
 */

void _square_tiling_integration( const char * _form_input, const char * _domain_input, const char *_output )
{
	Topology::CSTIMesh holo_mesh;
	Topology::CSTIMesh fund_mesh;

	holo_mesh.read_m( _form_input   );
	fund_mesh.read_m( _domain_input );

	Topology::CSquareTilingIntegration<Topology::CSTIMesh> integrator( & holo_mesh, & fund_mesh );
	integrator._integrate();

	fund_mesh.write_m( _output );
}

/*!	Linear combination of holomorphic 1-forms
*
*/

void _linear_combination( int argc, char * argv[] )
{

		Topology::CIMesh holo_mesh;
		double lambda;

		holo_mesh.read_m( argv[3] );
		lambda = atof( argv[4] );
	
		Topology::CLinearCombination<Topology::CIMesh> LC( &holo_mesh, lambda );

		for(int k = 5; k < argc; k ++ )
		{
			Topology::CIMesh form_mesh;
			double coefficient;
			
			form_mesh.read_m( argv[k] );
			k ++;
			coefficient = atof( argv[k] );
			LC._LinearCombine( &form_mesh, coefficient);
		}

		holo_mesh.write_m( argv[2] );

}

/*!	Linear combination of 1-forms
*
*/

void _linear_combination_one_form( int argc, char * argv[] )
{
		Holomorphy::CHCFMesh holo_mesh;
		double lambda;

		holo_mesh.read_m( argv[3] );
		lambda = atof( argv[4] );
	
		Topology::CLinearCombinationOneForm<Holomorphy::CHCFMesh> LC( &holo_mesh, lambda );

		for(int k = 5; k < argc; k ++ )
		{
			Holomorphy::CHCFMesh form_mesh;
			double coefficient;
			
			form_mesh.read_m( argv[k] );
			k ++;
			coefficient = atof( argv[k] );
			LC._LinearCombine( &form_mesh, coefficient);
		}

		holo_mesh.write_m( argv[2] );

}

/*!	compute the cut graph
 *
 */

void _cut_graph( const char * _input, const char * _output )
{
	Topology::CutGraphMesh spm;
	spm.read_m( _input );

	Topology::CCutGraph<Topology::CutGraphMesh> sp( & spm );
	sp._cut_locus();
	spm.write_m( _output );
};

/*!	compute the shortest path connecting an inner boundary to the exterior boundary
 *
 */
void _cut_domain( const char * _domain_mesh, const char * _mesh_with_cut )
{
	CSPMesh spm;
	spm.read_m( _domain_mesh );

	CShortestPath sp( & spm );
	sp._cut( _mesh_with_cut );
}


/*! Slice the mesh open along the sharp edges
 *
 */
void _slice( const char * _closed_mesh, const char * _open_mesh )
{
	Topology::CSMesh mesh;
	mesh.read_m( _closed_mesh );

	Topology::CWMesh wmesh( & mesh );
	wmesh.Slice();

	wmesh.wmesh()->write_m( _open_mesh );
};

/*!	compute the homology group basis
 *
 */
void _homology( const char * _input, const char * _output )
{
	Topology::CutGraphMesh spm;
	spm.read_m( _input );

	Topology::CHomology<Topology::CutGraphMesh> homology( & spm );
	homology._calculate_base();
	homology._output( _output );
};

/*!	compute the closed one form basis
 *
 */

void _cohomology_one_form( const char * _closed_mesh, const char * _open_mesh, const char * _closed_1_form, const char * _integrated )
{
	Topology::CCHMesh cmesh;
	cmesh.read_m( _closed_mesh );

	Topology::CCHMesh omesh;
	omesh.read_m( _open_mesh );

	Topology::CCohomology<Topology::CCHMesh> closed_form( & cmesh, &omesh );
	closed_form.calculate_closed_form();

	cmesh.write_m( _closed_1_form );
	omesh.write_m( _integrated );
}

/*!	compute the closed one form basis for multiply connected domain
 *
 */
void _cohomology_one_form_domain( const char * _input_mesh, const char * _sliced_mesh, const char * _closed_1_form, const char * _integrated_mesh )
{
	Topology::CCHMesh cmesh;
	cmesh.read_m( _input_mesh );

	Topology::CCHMesh omesh;
	omesh.read_m( _sliced_mesh );

	Topology::CDomainCohomology<Topology::CCHMesh> closed_form( & cmesh, &omesh );
	closed_form.calculate_closed_form();

	cmesh.write_m( _closed_1_form );
	omesh.write_m( _integrated_mesh );
}

/*! Puncture a hole in the center of the mesh
 *
 */
void _puncture( const char * _input_mesh, const char * _mesh_with_hole )
{
	Topology::CPMesh cmesh;
	cmesh.read_m( _input_mesh );

	Topology::CPuncture<Topology::CPMesh> cut( & cmesh );
	cut._puncture();

	cmesh.write_m( _mesh_with_hole );
}

/*! Fill the hole in the center
 *
 */
void _fill_puncture( const char * _mesh_with_hole, const char * _filled_mesh )
{

	Topology::CPMesh cmesh;
	cmesh.read_m( _mesh_with_hole );

	Topology::CPuncture<Topology::CPMesh> cut( & cmesh );
	cut._fill_hole();

	cmesh.write_m( _filled_mesh );
}


/*! Koebe's Iteration Method: Fill the central hole
 */
void Koebe_fill_centeral_hole( const char * _mesh_with_hole, const char * _filled_mesh )
{
	Topology::CSTMesh mesh;
	mesh.read_m( _mesh_with_hole );
	
	Topology::CStitch stitch( &mesh );
	stitch._fill_the_biggest_hole( _filled_mesh );
}

/*! Koebe's iteration Method: Remove faces with a given segment ID
 */
void Koebe_remove_segment( const char * _input_mesh, const char * segment_id, const char * _output_mesh )
{

	Topology::CSTMesh mesh;
	mesh.read_m( _input_mesh );
	
	Topology::CStitch stitch( &mesh );
	int sid = atoi( segment_id );
	stitch._remove_segment(sid,  _output_mesh );
}

/*! Koebe's iteration method: Copy UV, from UV_Mesh to Position_Mesh
 */
void Koebe_copy_uv( const char * _uv_mesh, const char * _pos_mesh, const char * _output_mesh )
{
	Topology::CSTMesh uv_mesh;
	uv_mesh.read_m( _uv_mesh );

	Topology::CSTMesh position_mesh;
	position_mesh.read_m( _pos_mesh );
	
	Topology::CStitch stitch( &uv_mesh );
	stitch._copy_uv( &position_mesh );
	position_mesh.write_m( _output_mesh );

}

/*! Compute the double covering mesh
 */
void _double_covering( const char * _open_mesh, const char * _doubled_mesh )
{
		Topology::CDCMesh dmesh;
		dmesh.read_m( _open_mesh );
		dmesh.DoubleCovering();
		dmesh.write_m( _doubled_mesh );
}


/*! Compute the shortest basis using Erickson's method
 */
void _shortest_homology_basis( const char * _input, const char * _vertex_id, const char * _basis_mesh )
{

	CDKMesh mesh;
	mesh.read_m( _input );
	CDijkstra ST( &mesh, &mesh );
	CDijkstraVertex * pV = mesh.idVertex( atoi(_vertex_id ) );
	ST._shortest_basis( pV );	
	mesh.write_m( _basis_mesh );
}


