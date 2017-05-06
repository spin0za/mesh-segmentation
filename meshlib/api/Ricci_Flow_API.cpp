#include "Ricci_Flow_API.h"

using namespace MeshLib::RicciFlow;

unsigned int CRicciFlowVertex::traits = 0;


/***************************************************************************************************************************

	For topological Torus


****************************************************************************************************************************/

//-genus_one_tangential_ricci_flow gnome.m gnome.open.m gnome.uv.m
void _genus_one_tangential_ricci_flow( const char * _closed_mesh, const char * _fundamental_domain, const char * _mesh_with_uv )
{

	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CRFMesh mesh;
	mesh.read_m( _closed_mesh );
	CGenusOneTangentialRicciFlow<CRFMesh> mapper(&mesh);

	mapper._calculate_metric();

	CRFMesh omesh;
	omesh.read_m( _fundamental_domain );
	
	for( CRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		CRicciFlowEdge * pE = *eiter;
		CRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		CRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
		CRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
		CRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
		CRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	CRFEmbed embed( &omesh );
	embed._embed();
	omesh.write_m( _mesh_with_uv );

}

//-ricci_flow kitten.m kitten.open.m kitten.uv.m
void _tangential_ricci_flow( const char * _closed_mesh, const char * _fundamental_domain, const char * _mesh_with_uv )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CRFMesh mesh;
	mesh.read_m( _closed_mesh );
	CTangentialRicciFlow<CRFMesh> mapper(&mesh);

	mapper._calculate_metric();

	CRFMesh omesh;
	omesh.read_m( _fundamental_domain );
	
	for( CRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		CRicciFlowEdge * pE = *eiter;
		CRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		CRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
		CRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
		CRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
		CRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	CRFEmbed embed( &omesh );
	embed._embed();
	omesh.write_m( _mesh_with_uv );
}


//-ricci_flow sophie.remesh.m sophie.uv.m
void _inverse_ricci_flow( const char * _closed_mesh, const char * _fundamental_domain, const char * _mesh_with_uv )
{
		CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

		CRFMesh mesh;
		mesh.read_m( _closed_mesh );

		CInversiveDistanceRicciFlow<CRFMesh> mapper(&mesh);
		mapper._calculate_metric();

		CRFMesh omesh;
		omesh.read_m( _fundamental_domain );
		
		for( CRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
		{
			CRicciFlowEdge * pE = *eiter;
			CRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
			CRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
			CRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
			CRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
			CRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
			pE->length() = pWe->length();
		}

		CRFEmbed embed( &omesh );
		embed._embed();
		omesh.write_m( _mesh_with_uv );
}

//-yamabe_flow kitten.m kitten.open.m kitten.uv.m
void _Yamabe_flow( const char * _closed_mesh, const char * _fundamental_domain, const char * _mesh_with_uv )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CRFMesh mesh;
	mesh.read_m( _closed_mesh );
	CYamabeFlow<CRFMesh> mapper(&mesh);

	mapper._calculate_metric();

	CRFMesh omesh;
	omesh.read_m( _fundamental_domain );
	
	for( CRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		CRicciFlowEdge * pE = *eiter;
		CRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		CRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
		CRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
		CRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
		CRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	CRFEmbed embed( &omesh );
	embed._embed();
	omesh.write_m( _mesh_with_uv );

}

/*
 * -virtual_inversive_ricci_flow kitten.m kitten.open.m kitten.uv.m
 */
void _virtual_inversive_ricci_flow( const char * input_mesh, const char * fundamental_domain, const char * uv_mesh )
{
		CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

		CRFMesh mesh;
		mesh.read_m( input_mesh );

		CVirtualInversiveDistanceRicciFlow<CRFMesh> mapper(&mesh);
		mapper._calculate_metric( 0.25 );

		CRFMesh omesh;
		omesh.read_m( fundamental_domain );
		
		for( CRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
		{
			CRicciFlowEdge * pE = *eiter;
			CRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
			CRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
			CRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
			CRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
			CRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
			pE->length() = pWe->length();
		}

		CRFEmbed embed( &omesh );
		embed._embed();
		omesh.write_m( uv_mesh );
}

/******************************************************************************************************************************
*
*	Extremal Length
*
*******************************************************************************************************************************/

//-yamabe_extremal_length sophie.remesh.m sophie.uv.m
void _Yamabe_extremal_length( const char * _input_mesh, const char * _mesh_with_uv )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CRFMesh mesh;
	mesh.read_m( _input_mesh );

	CYamabeExtremalLength<CRFMesh> mapper(&mesh);
	mapper._calculate_metric();

	CRFEmbed embed( &mesh );
	embed._embed();
	mesh.write_m( _mesh_with_uv );
}

//-idrf_extremal_length sophie.remesh.m sophie.uv.m
void _idrf_extremal_length( const char * _input_mesh, const char * _mesh_with_uv )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CRFMesh mesh;
	mesh.read_m( _input_mesh );

	CInversiveDistanceRicciFlowExtremalLength<CRFMesh> mapper(&mesh);
	mapper._calculate_metric();

	CRFEmbed embed( &mesh );
	embed._embed();
	mesh.write_m( _mesh_with_uv );
}

//-tangent_ricci_extremal_length sophie.remesh.m sophie.uv.m
void _tangent_ricci_extremal_length( const char * _input_mesh, const char * _mesh_with_uv )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CRFMesh mesh;
	mesh.read_m( _input_mesh );

	CTangentialRicciFlowExtremalLength<CRFMesh> mapper(&mesh);
	mapper._calculate_metric();

	CRFEmbed embed( &mesh );
	embed._embed();
	mesh.write_m( _mesh_with_uv );
}
//-virtual_extremal_length sophie.m sophie.uv.m
void _virtual_inversive_distance_extremal_length( const char * input_mesh, const char * uv_mesh )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CRFMesh mesh;
	mesh.read_m( input_mesh );

	CVirtualInversiveDistanceRicciFlowExtremalLength<CRFMesh> mapper(&mesh);
	mapper._calculate_metric( 0.25 );

	CRFEmbed embed( &mesh );
	embed._embed();
	mesh.write_m( uv_mesh );
}

/******************************************************************************************************************************
*
*	Riemann Mapping
*
*******************************************************************************************************************************/


//-ricci_riemann_map sophie.m sophie.open.m sophie.uv.m
void _tangent_ricci_riemann_map( const char * _closed_mesh, const char * _fundamental_domain, const char * _mesh_with_uv )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CRFMesh mesh;
	mesh.read_m( _closed_mesh );
	//CTangentialRicciFlow<CRicciFlowVertex,CRicciFlowEdge,CRicciFlowFace,CRicciFlowHalfEdge> mapper(&mesh);
	
	CTangentialRicciFlowRiemannMapping<CRFMesh> mapper(&mesh);		
	mapper._calculate_metric();

	CRFMesh omesh;
	omesh.read_m( _fundamental_domain );
	
	for( CRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		CRicciFlowEdge * pE = *eiter;
		CRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		CRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
		CRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
		CRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
		CRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	CRFEmbed embed( &omesh );
	embed._embed();

	CCurvatureFlowExponentialMap<CRFMesh> exp_mapper(&mesh, &omesh);

	mesh.write_m( _mesh_with_uv );
}


//-idrf_riemann_map sophie.m sophie.open.m sophie.uv.m
void _idrf_riemann_map(  const char * _closed_mesh, const char * _fundamental_domain, const char * _mesh_with_uv )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CRFMesh mesh;
	mesh.read_m( _closed_mesh );
	
	CInversiveDistanceRicciFlowRiemannMapping<CRFMesh> mapper(&mesh);		
	mapper._calculate_metric();

	CRFMesh omesh;
	omesh.read_m( _fundamental_domain );
	
	for( CRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		CRicciFlowEdge * pE = *eiter;
		CRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		CRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
		CRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
		CRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
		CRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	CRFEmbed embed( &omesh );
	embed._embed();
	CCurvatureFlowExponentialMap<CRFMesh> exp_mapper(&mesh, &omesh);
	mesh.write_m( _mesh_with_uv );
}


//-yamabe_riemann_map sophie.m sophie.open.m sophie.uv.m
void _yamabe_riemann_map( const char * _closed_mesh, const char * _fundamental_domain, const char * _mesh_with_uv )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CRFMesh mesh;
	mesh.read_m( _closed_mesh );
	
	CInversiveDistanceRicciFlowRiemannMapping<CRFMesh> mapper(&mesh);		
	mapper._calculate_metric();

	CRFMesh omesh;
	omesh.read_m( _fundamental_domain );
	
	for( CRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		CRicciFlowEdge * pE = *eiter;
		CRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		CRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
		CRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
		CRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
		CRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	CRFEmbed embed( &omesh );
	embed._embed();
	CCurvatureFlowExponentialMap<CRFMesh> exp_mapper(&mesh, &omesh);
	mesh.write_m( _mesh_with_uv );

}


/*!
 * -virtual_riemann_mapping sophie.m sophie.open.m sophie.uv.m
 */
void _virtual_distance_riemann_mapping( const char * input_mesh, const char * open_mesh, const char * uv_mesh )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CRFMesh mesh;
	mesh.read_m( input_mesh );

	CVirtualInversiveDistanceRicciFlowRiemannMapping<CRFMesh> mapper(&mesh);
	mapper._calculate_metric(0.25);

	CRFMesh omesh;
	omesh.read_m( open_mesh );
	
	for( CRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		CRicciFlowEdge * pE = *eiter;
		CRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		CRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
		CRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
		CRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
		CRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	CRFEmbed embed( &omesh );
	embed._embed();

	CCurvatureFlowExponentialMap<CRFMesh> exp_mapper(&mesh, &omesh);

	mesh.write_m( uv_mesh );

}


/******************************************************************************************************************************
*
*	Poly Annulus
*
*******************************************************************************************************************************/

//-ricci_poly_annulus alex_hole.m alex_hole.open.m alex_hole.uv.m
void _ricci_poly_annulus( const char * _input_mesh, const char * _fundamental_domain, const char * _mesh_with_uv )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CRFMesh mesh;
	mesh.read_m( _input_mesh );
	
	CTangentialRicciFlowPolyAnnulus<CRFMesh> mapper(&mesh);		
	mapper._calculate_metric();

	CRFMesh omesh;
	omesh.read_m( _fundamental_domain );
	
	for( CRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		CRicciFlowEdge * pE = *eiter;
		CRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		CRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
		CRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
		CRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
		CRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	CRFEmbed embed( &omesh );
	embed._embed();
	//omesh.write_m( argv[4] );

	CCurvatureFlowExponentialMap<CRFMesh> exp_mapper(&mesh, &omesh);
	mesh.write_m( _mesh_with_uv );
}

//-idrf_poly_annulus sophie.m sophie.open.m sophie.uv.m
void _idrf_poly_annulus( const char * _input_mesh, const char * _fundamental_domain, const char * _mesh_with_uv )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CRFMesh mesh;
	mesh.read_m( _input_mesh );
	
	CInversiveDistanceRicciFlowPolyAnnulus<CRFMesh> mapper(&mesh);		
	mapper._calculate_metric();

	CRFMesh omesh;
	omesh.read_m( _fundamental_domain );
	
	for( CRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		CRicciFlowEdge * pE = *eiter;
		CRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		CRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
		CRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
		CRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
		CRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	CRFEmbed embed( &omesh );
	embed._embed();
	//omesh.write_m( argv[4] );

	CCurvatureFlowExponentialMap<CRFMesh> exp_mapper(&mesh, &omesh);

	mesh.write_m( _mesh_with_uv );
}


//-yamabe_poly_annulus sophie.m sophie.open.m sophie.uv.m
void _yamabe_poly_annulus( const char * _input_mesh, const char * _fundamental_domain, const char * _mesh_with_uv )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CRFMesh mesh;
	mesh.read_m( _input_mesh );
	
	CYamabeFlowPolyAnnulus<CRFMesh> mapper(&mesh);		
	mapper._calculate_metric();

	CRFMesh omesh;
	omesh.read_m( _fundamental_domain );
	
	for( CRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		CRicciFlowEdge * pE = *eiter;
		CRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		CRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
		CRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
		CRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
		CRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	CRFEmbed embed( &omesh );
	embed._embed();
	//omesh.write_m( argv[4] );

	CCurvatureFlowExponentialMap<CRFMesh> exp_mapper(&mesh, &omesh);

	mesh.write_m( _mesh_with_uv );
}

/*!
 * -virtual_poly_annulus alex_hole.m alex_hole.open.m alex_hole.uv.m
 */
void _virtual_distance_poly_annulus( const char * input_mesh, const char * open_mesh, const char * uv_mesh )
{
	//needs to save vertex uv trait

	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CRFMesh mesh;
	mesh.read_m( input_mesh );
	
	CVirtualInversiveDistanceRicciFlowPolyAnnulus<CRFMesh> mapper(&mesh);		
	mapper._calculate_metric();

	CRFMesh omesh;
	omesh.read_m( open_mesh );
	
	for( CRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		CRicciFlowEdge * pE = *eiter;
		CRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		CRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
		CRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
		CRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
		CRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	CRFEmbed embed( &omesh );
	embed._embed();
	//omesh.write_m( argv[4] );

	CCurvatureFlowExponentialMap<CRFMesh> exp_mapper(&mesh, &omesh);
	mesh.write_m( uv_mesh );
}


/******************************************************************************************************************************
*
*	Hyperbolic Ricci Flow
*
*******************************************************************************************************************************/

/*!
 * -hyper_ricci_flow eight.m eight.open.m eight.uv.m
 */
void _hyperbolic_tangential_ricci_flow( const char * closed_mesh, const char * fundamental_domain, const char * embedded_mesh )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CHRFMesh mesh;
	mesh.read_m( closed_mesh );
	CTangentialHyperbolicRicciFlow<CHRFMesh> mapper(&mesh);

	mapper._calculate_metric();

	CRFMesh omesh;
	omesh.read_m( fundamental_domain );
	
	for( CRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		CRicciFlowEdge * pE = *eiter;
		CRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		CRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
		CRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
		CRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
		CRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	CHRFEmbed embed( &omesh );
	embed._embed();
	omesh.write_m( embedded_mesh );
}


//-hyper_ricci_flow eight.m eight.open.m eight.uv.m
void _hyperbolic_inverse_ricci_flow( const char * closed_mesh, const char * fundamental_domain, const char * embedded_mesh )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CHRFMesh mesh;
	mesh.read_m( closed_mesh );
	CInversiveDistanceHyperbolicRicciFlow<CHRFMesh> mapper(&mesh);

	mapper._calculate_metric();

	CRFMesh omesh;
	omesh.read_m( fundamental_domain );
	
	for( CRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		CRicciFlowEdge * pE = *eiter;
		CRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		CRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
		CRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
		CRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
		CRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	CHRFEmbed embed( &omesh );
	embed._embed();
	omesh.write_m( embedded_mesh );

}

/*!
 * -hyper_yamabe_flow eight.m eight.open.m eight.uv.m
 */
void _hyperbolic_yamabe_flow( const char * closed_mesh, const char * open_mesh, const char * uv_mesh )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CHRFMesh mesh;
	mesh.read_m( closed_mesh );
	CHyperbolicYamabeFlow<CHRFMesh> mapper(&mesh);

	mapper._calculate_metric(1e-6, 0.05);

	CRFMesh omesh;
	omesh.read_m( open_mesh );
	
	for( CRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		CRicciFlowEdge * pE = *eiter;
		CRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		CRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
		CRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
		CRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
		CRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	CHRFEmbed embed( &omesh );
	embed._embed();
	omesh.write_m( uv_mesh );

}

/******************************************************************************************************************************
*
*	Tangential Hyperbolic RIcci Flow for Riemann Mapping
*
*******************************************************************************************************************************/


/*!
 * set boundary vertex radius to go to infinity, -hyper_ricci_flow_riemann_map face.m 2.0 face.uv.m
 */
void _hyperbolic_ricci_flow_riemann_mapping( const char * input_mesh, const char * boundary_radius, const char * uv_mesh, const char * output_eps )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CHRFMesh mesh;
	mesh.read_m( input_mesh );
	CTangentialHyperbolicRicciFlowRiemannMap<CHRFMesh> mapper(&mesh, atof( boundary_radius ));

	mapper._calculate_metric();

	CRFMesh omesh;
	omesh.read_m( input_mesh );
	
	for( CRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		CRicciFlowEdge * pE = *eiter;
		CRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		CRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
		CRicciFlowVertex * pW0 = mesh.idVertex( pV0->id() );
		CRicciFlowVertex * pW1 = mesh.idVertex( pV1->id() );
		CRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	for( CRFMesh::MeshVertexIterator viter( &omesh ); !viter.end(); viter ++ )
	{
		CRicciFlowVertex * pV = *viter;
		CRicciFlowVertex * pW = mesh.idVertex( pV->id() );
		pV->u() = pW->u();
	}

	CHRFEmbed embed( &omesh );
	embed._embed();
	omesh.write_m( uv_mesh );

	CPs ps( &omesh );

	ps.hyperbolic_print( output_eps );

}

/******************************************************************************************************************************
*
*	Quasi-Conformal Ricci Flow Riemann Mapping
*
*******************************************************************************************************************************/
/*
 * -qc_yamabe_flow_riemann_map sophie.m sophie.open.m sophie.uv.m
 */
void _qc_yamabe_riemann_mapping( const char * input_mesh, const char * fundamental_domain, const char * uv_mesh )
{

	CQCRFMesh mesh;
	mesh.read_m( input_mesh );
	CQCYamabeFlowRiemannMapping<CQCRFMesh> mapper(&mesh);
	mapper._calculate_metric();

	CQCRFMesh omesh;
	omesh.read_m( fundamental_domain );
	
	for( CQCRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		CQCRicciFlowEdge * pE = *eiter;
		CQCRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		CQCRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
		CQCRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
		CQCRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
		CQCRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	
	CEuclideanEmbed<CQCRFMesh> embed( &omesh );
	embed._embed();

	CCurvatureFlowExponentialMap<CQCRFMesh> exp_mapper( &mesh, &omesh);

	mesh.write_m( uv_mesh );

}

/*
 * -qc_inversive_distance_riemann_map sophie.m sophie.open.m sophie.uv.m
 */
void _qc_inversive_distance_riemann_mapping( const char * input_mesh, const char * fundamental_domain, const char * uv_mesh )
{

	CQCRFMesh mesh;
	mesh.read_m( input_mesh );
	CQCInversiveDistanceRicciFlow<CQCRFMesh> mapper(&mesh);
	mapper._calculate_metric();

	CQCRFMesh omesh;
	omesh.read_m( fundamental_domain );
	
	for( CQCRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		CQCRicciFlowEdge * pE = *eiter;
		CQCRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		CQCRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
		CQCRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
		CQCRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
		CQCRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	
	CEuclideanEmbed<CQCRFMesh> embed( &omesh );
	embed._embed();

	CCurvatureFlowExponentialMap<CQCRFMesh> exp_mapper( &mesh, &omesh);

	mesh.write_m( uv_mesh );

}

/*
 * -qc_virtual_inversive_distance_riemann_map sophie.m sophie.open.m sophie.uv.m
 */
//if( strcmp( argv[1] , "-qc_ricci_riemann_map") == 0 && argc == 5)
void _qc_virtual_inversive_distance_riemann_mapping( const char * input_mesh, const char * fundamental_domain, const char * uv_mesh )
{

	CQCRFMesh mesh;
	mesh.read_m( input_mesh );
	CQCVirtualInversiveDistanceRicciFlow<CQCRFMesh> mapper(&mesh);
	mapper._calculate_metric( 0.25 );

	CQCRFMesh omesh;
	omesh.read_m( fundamental_domain );
	
	for( CQCRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		CQCRicciFlowEdge * pE = *eiter;
		CQCRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		CQCRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
		CQCRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
		CQCRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
		CQCRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	
	CEuclideanEmbed<CQCRFMesh> embed( &omesh );
	embed._embed();

	CCurvatureFlowExponentialMap<CQCRFMesh> exp_mapper( &mesh, &omesh);

	mesh.write_m( uv_mesh );

}


//-dynamic_yamabe_flow kitten.m kitten.uv.m
void _dynamic_Yamabe_flow( const char * _closed_mesh, const char * _mesh_with_uv )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CDRFMesh mesh;
	mesh.read_m( _closed_mesh );

/* for Debug purpose

	for( CDRFMesh::MeshEdgeIterator eiter( &mesh ); !eiter.end(); eiter ++ )
	{
		CRicciFlowEdge * pE = *eiter;
		pE->length() = mesh.edgeLength( pE );
	}

	mesh.write_g("initial.g");
	mesh._preserve_Delaunay();
*/

	CDynamicYamabeFlow<CDRFMesh> mapper(&mesh);

	mapper._calculate_metric();
/*
	CRFMesh omesh;
	omesh.read_m( _fundamental_domain );
	
	for( CRFMesh::MeshEdgeIterator eiter( &omesh ); !eiter.end(); eiter ++ )
	{
		CRicciFlowEdge * pE = *eiter;
		CRicciFlowVertex * pV0 = omesh.edgeVertex1( pE );
		CRicciFlowVertex * pV1 = omesh.edgeVertex2( pE );
		CRicciFlowVertex * pW0 = mesh.idVertex( pV0->father() );
		CRicciFlowVertex * pW1 = mesh.idVertex( pV1->father() );
		CRicciFlowEdge * pWe = mesh.vertexEdge( pW0, pW1 );
		pE->length() = pWe->length();
	}

	CRFEmbed embed( &omesh );
	embed._embed();
	omesh.write_m( _mesh_with_uv );
*/
}


