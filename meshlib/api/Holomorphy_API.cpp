#include "Holomorphy_API.h"
#include "Ricci_Flow_API.h"

using namespace MeshLib;


//compute harmonic map, between a topological disk to a disk
void _harmonic_map( const char * _input, const char * _output )
{
	Holomorphy::CHMMesh mesh;
	mesh.read_m( _input );

	Holomorphy::CDiskHarmonicMapper<Holomorphy::CHMMesh> mapper( & mesh );
	mapper._map();
	mesh.write_m( _output );
};


//compute spherical harmonic map from the topological sphere to the unit sphere

void _spherical_harmonic_map( const char * _input, const char * _output )
{
	Holomorphy::CSHMMesh mesh;
	mesh.read_m( _input );

	Holomorphy::CSphericalHarmonicMap<Holomorphy::CSHMMesh> mapper( & mesh );
	mapper.map();

	for( Holomorphy::CSHMMesh::MeshVertexIterator viter( &mesh ); !viter.end(); ++ viter )
	{
		Holomorphy::CSHMMesh::CVertex * v = *viter;
		v->point() = v->u();
	}

	mesh.write_m( _output );
}

/*!	compute the harmonic form for extremal length, need to compare to holomorphic form, and integrate
 *
 */

void _extremal_length( const char * _input, const char * _output )
{
	Holomorphy::CELMesh cmesh;
	cmesh.read_m( _input );


	Holomorphy::CExtremalLength<Holomorphy::CELMesh> EL( & cmesh );
	EL.calculate_harmonic_exact_form( _output );
}

/*!	compute the harmonic form for extremal length, using combinatorial Laplacian
 *
 */

void _combinatorial_extremal_length( const char * _input, const char * _output )
{
	Holomorphy::CELMesh cmesh;
	cmesh.read_m( _input );


	Holomorphy::CCombinatorialExtremalLength<Holomorphy::CELMesh> EL( & cmesh );
	EL.calculate_harmonic_exact_form( _output );
}

/*!	Compute the slit maps
*
*/
void _slit_map( int argc, char * argv[] )
{

	std::vector<Holomorphy::CSMMesh*> meshes;

	printf("meshes size %d\n", meshes.size() );

	for( int i = 2; i < argc-3; i ++ )
	{
		Holomorphy::CSMMesh * pMesh = new Holomorphy::CSMMesh;
		assert( pMesh );
		pMesh->read_m( argv[i] );
		meshes.push_back( pMesh );
	}

	Holomorphy::CSlitMap<Holomorphy::CSMMesh> map( meshes );

	int c1 = atoi( argv[argc-3] );
	int c2 = atoi( argv[argc-2] );

	map._slit_map( c1,c2 );

	meshes[0]->write_m( argv[argc-1] );

	for( size_t i = 0; i < meshes.size(); i ++ )	
	{
		delete meshes[i];
	}
}


/*! compute the basis for all the exact harmonic forms on a surface
 *
 */
void _exact_form( const char * _input, const char * _output )
{
	Holomorphy::CHarmonicMesh hm;
	hm.read_m( _input );

	Holomorphy::CHarmonicExactForm<Holomorphy::CHarmonicMesh> exact_form( & hm );
	exact_form.calculate_harmonic_exact_form( _output );
}

/*!	compute the basis of all holomorphic 1-forms on a surface
 *
 */
void _holomorphic_form( int argc,  char * argv[] )
{
	std::list<Holomorphy::CHoloFormMesh*> meshes;


	for( int i = 2; i < argc; i ++ )
	{
		Holomorphy::CHoloFormMesh * pMesh = new Holomorphy::CHoloFormMesh;
		assert( pMesh );
		pMesh->read_m( argv[i] );
		meshes.push_back( pMesh );
	}

	Holomorphy::CHolomorphicForm<Holomorphy::CHoloFormMesh> form( meshes );
	form.conjugate();

	int id = 0;
	for( std::list<Holomorphy::CHoloFormMesh*>::iterator miter = meshes.begin(); miter != meshes.end(); miter++)
	{
		Holomorphy::CHoloFormMesh * pM = *miter;
		std::stringstream iss;
		iss << "holomorphic_form_" << id++ << ".m";
		pM->write_m( iss.str().c_str() );
	}

	for( std::list<Holomorphy::CHoloFormMesh*>::iterator miter = meshes.begin(); miter != meshes.end(); miter++)
	{
		Holomorphy::CHoloFormMesh * pM = *miter;
		delete pM;
	}
}

/*! Compute the exponential map
 *
 */
void _polar_map( const char * _closed_mesh, const char * _open_mesh, const char * _output )
{
	Holomorphy::CPMMesh mesh;
	mesh.read_m( _closed_mesh );

	Holomorphy::CPMMesh open_mesh;
	open_mesh.read_m( _open_mesh );

	Holomorphy::CPolarMap<Holomorphy::CPMMesh> map( &mesh, &open_mesh );
	map._exponential_map();

	mesh.write_m( _output );
}

/*!	diffuse the closed one form basis to a harmonic one form
 *
 */
void _diffuse( const char * _closed_1_form, const char * _fundamental_domain, const char * _harmonic_1_form, const char * _harmonic_form_integrated )
{
	Holomorphy::CHCFMesh cmesh;
	cmesh.read_m( _closed_1_form );

	Holomorphy::CHCFMesh omesh;
	omesh.read_m( _fundamental_domain );

	Holomorphy::CHeatFlow<Holomorphy::CHCFMesh> flow( & cmesh, &omesh );
	flow.calculate_harmonic_form();

	cmesh.write_m( _harmonic_1_form );
	omesh.write_m( _harmonic_form_integrated );
}

/*!	diffuse the closed one form basis to a harmonic one form for square tiling purpose
 *
 */
void _combinatorial_diffuse( const char * _closed_1_form, const char * _fundamental_domain, const char * _harmonic_1_form, const char * _harmonic_form_integrated )
{
	Holomorphy::CHCFMesh cmesh;
	cmesh.read_m( _closed_1_form );

	Holomorphy::CHCFMesh omesh;
	omesh.read_m( _fundamental_domain );

	Holomorphy::CCombinatorialHeatFlow<Holomorphy::CHCFMesh> flow( & cmesh, &omesh );
	flow.calculate_harmonic_form();

	cmesh.write_m( _harmonic_1_form );
	omesh.write_m( _harmonic_form_integrated );
}


/******************************************************************************************************************************************
*
*
*	Quasi-Conformal Geometric Methods
*
*
*******************************************************************************************************************************************/

/*! compute the basis for all the exact harmonic forms on a surface
 *
 */
void _qc_exact_form( const char * _input, const char * _output )
{
	Holomorphy::CHarmonicMesh hm;
	hm.read_m( _input );
		
	_read_vertex_z<Holomorphy::CHarmonicMesh, Holomorphy::CHVertex, Holomorphy::CHEdge, CFace, Holomorphy::CHHalfEdge>( &hm );
	_read_vertex_mu<Holomorphy::CHarmonicMesh, Holomorphy::CHVertex, Holomorphy::CHEdge, CFace, Holomorphy::CHHalfEdge>( &hm );

	Holomorphy::CQCHarmonicExactForm<Holomorphy::CHarmonicMesh> exact_form( & hm );
	exact_form.calculate_harmonic_exact_form( _output );
}


/*!	diffuse the closed one form basis to a quasi-conformal harmonic one form
 *
 */

void _qc_diffuse( const char * _closed_1_form, const char * _fundamental_domain, const char * _harmonic_1_form, const char * _integrated )
{
		Holomorphy::CHCFMesh cmesh;
		cmesh.read_m( _closed_1_form );

		_read_vertex_z<Holomorphy::CHCFMesh, Holomorphy::CHCFVertex, Holomorphy::CHCFEdge, CFace, Holomorphy::CHCFHalfEdge>( &cmesh );
		_read_vertex_mu<Holomorphy::CHCFMesh, Holomorphy::CHCFVertex, Holomorphy::CHCFEdge, CFace, Holomorphy::CHCFHalfEdge>( &cmesh );

		Holomorphy::CHCFMesh omesh;
		omesh.read_m( _fundamental_domain );

		_read_vertex_z<Holomorphy::CHCFMesh, Holomorphy::CHCFVertex, Holomorphy::CHCFEdge, CFace, Holomorphy::CHCFHalfEdge>( &omesh );
		_read_vertex_mu<Holomorphy::CHCFMesh, Holomorphy::CHCFVertex, Holomorphy::CHCFEdge, CFace, Holomorphy::CHCFHalfEdge>( &omesh );

		Holomorphy::CQCHeatFlow<Holomorphy::CHCFMesh> flow( & cmesh, &omesh );
		
		flow.calculate_harmonic_form();
		cmesh.write_m( _harmonic_1_form );
		omesh.write_m( _integrated );
}

/*!	compute the basis of all holomorphic 1-forms on a surface
 *
 */
void _qc_holomorphic_form( int argc, char * argv[] )
{
	std::list<Holomorphy::CHoloFormMesh*> meshes;


	for( int i = 2; i < argc; i ++ )
	{
		Holomorphy::CHoloFormMesh * pMesh = new Holomorphy::CHoloFormMesh;
		assert( pMesh );
		pMesh->read_m( argv[i] );
		meshes.push_back( pMesh );

		_read_vertex_z<Holomorphy::CHoloFormMesh,Holomorphy::CHoloFormVertex, Holomorphy::CHoloFormEdge, Holomorphy::CHoloFormFace, Holomorphy::CHoloFormHalfEdge>( pMesh );
		_read_vertex_mu<Holomorphy::CHoloFormMesh,Holomorphy::CHoloFormVertex, Holomorphy::CHoloFormEdge, Holomorphy::CHoloFormFace, Holomorphy::CHoloFormHalfEdge>( pMesh );

	}

	Holomorphy::CQCHolomorphicForm<Holomorphy::CHoloFormMesh> form( meshes );
	form.conjugate();

	int id = 0;
	for( std::list<Holomorphy::CHoloFormMesh*>::iterator miter = meshes.begin(); miter != meshes.end(); miter++)
	{
		Holomorphy::CHoloFormMesh * pM = *miter;
		std::stringstream iss;
		iss << "holomorphic_form_" << id++ << ".m";
		pM->write_m( iss.str().c_str() );
	}

	for( std::list<Holomorphy::CHoloFormMesh*>::iterator miter = meshes.begin(); miter != meshes.end(); miter++)
	{
		Holomorphy::CHoloFormMesh * pM = *miter;
		delete pM;
	}
}

/**********************************************************************************************************************************************
*
*	Conformal Mapping Based on Diagnal Ratio
*	
*
**********************************************************************************************************************************************/

/*!	compute the basis of all holomorphic 1-forms based on diagonal ratio
 *
 */

void _diagonal_ratio_holomorphic_form( int argc, char * argv[] )
{
	std::list<Holomorphy::CDRMesh*> meshes;

	for( int i = 2; i < argc; i ++ )
	{
		Holomorphy::CDRMesh * pMesh = new Holomorphy::CDRMesh;
		assert( pMesh );
		pMesh->read_m( argv[i] );
		meshes.push_back( pMesh );
	}

	Holomorphy::CDiagonalRatioHolomorphicForm<Holomorphy::CDRMesh> form( meshes );
	form.conjugate();

	int id = 0;
	for( std::list<Holomorphy::CDRMesh*>::iterator miter = meshes.begin(); miter != meshes.end(); miter++)
	{
		Holomorphy::CDRMesh * pM = *miter;
		std::stringstream iss;
		iss << "holomorphic_form_" << id++ << ".m";
		pM->write_m( iss.str().c_str() );
	}

	for( std::list<Holomorphy::CDRMesh*>::iterator miter = meshes.begin(); miter != meshes.end(); miter++)
	{
		Holomorphy::CDRMesh * pM = *miter;
		delete pM;
	}
}


/*!	compute the conformal mapping based on diagonal ratio
 *
 */

void _diagonal_ratio_free_boundary_mapping( const char * _input_mesh, const char * _output_mesh )
{
	Holomorphy::CDRMesh mesh;
	mesh.read_m( _input_mesh );

	Holomorphy::CDiagonalRatioMapping<Holomorphy::CDRMesh> map( &mesh );
	map.map();
	mesh.write_m( _output_mesh );
}



/*!	compute the basis of all holomorphic quadratic forms
 *
 */
void _holomorphic_quadratic_form( int argc,  char * argv[] )
{
	std::list<Holomorphy::CHolomorphicQuadraticForm<Holomorphy::CHoloQuadFormMesh>*> forms;
	std::list<Holomorphy::CHoloQuadFormMesh*> meshes;


	for( int i = 2; i < argc; i ++ )
	{
		Holomorphy::CHoloQuadFormMesh * pMesh = new Holomorphy::CHoloQuadFormMesh;
		assert( pMesh );
		pMesh->read_m( argv[i] );
		meshes.push_back( pMesh );
		
		Holomorphy::CHolomorphicQuadraticForm<Holomorphy::CHoloQuadFormMesh> * form = new Holomorphy::CHolomorphicQuadraticForm<Holomorphy::CHoloQuadFormMesh>( pMesh );
		assert( form );
		form->normalize();
		forms.push_back( form );
	}

	std::vector<double> coefficients;

	coefficients.push_back( 1.0 );
	coefficients.push_back( 0.5 );

	forms.front()->linear_combine( meshes, coefficients );
	forms.front()->form_to_metric();
	
	Holomorphy::CHoloQuadFormMesh *  qm = meshes.front();


	RicciFlow::CRFMesh omesh;
	omesh.read_m( argv[2] );
	
	for( RicciFlow::CRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		RicciFlow::CRicciFlowEdge * pE = *eiter;
		RicciFlow::CRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		RicciFlow::CRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );

		Holomorphy::CHoloQuadFormVertex * pW0 = qm->idVertex( pV0->id() );
		Holomorphy::CHoloQuadFormVertex * pW1 = qm->idVertex( pV1->id() );

		Holomorphy::CHoloQuadFormEdge   * pWe = qm->vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	RicciFlow::CRicciFlowVertex::traits = RicciFlow::CRicciFlowVertex::traits | TRAIT_UV;

	RicciFlow::CRFEmbed embed( &omesh );
	embed._embed();
	omesh.write_m( "output.m" );
}

