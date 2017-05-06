#include "Delaunay_API.h"

using namespace MeshLib;

bool CDelaunayRemeshVertex::m_with_rgb=false;
bool CQTDelaunayRemeshVertex::m_with_rgb=false;

/*! Planar Mesh Generation
*
*/

void _Rupert_Delaunay_Refinement( const char * poly_file, const char * mesh_file )
{
    CDTMesh dmesh;
	dmesh.ProcessPoly( poly_file, mesh_file );
}

/*! Surface Remeshing based on Delaunay Refinement Method
*
*/
void _Surface_Delaunay_Remesh( const char * _input_mesh, const char * _output_mesh )
{
    
		CSDRMesh dmesh;
		dmesh.read_m( _input_mesh );

		for( CSDRMesh::MeshVertexIterator viter( &dmesh ); !viter.end(); viter ++ )
		{
			CDelaunayRemeshVertex * pV = *viter;
			pV->uv() = pV->uv() * 2.0 - CPoint2(1,1);
		}
		dmesh._lambda();

		CSDRMesh tmesh;
		tmesh.Remesh( dmesh );

		CSDRMesh omesh;
		tmesh.Convert( omesh );
		omesh.write_m( _output_mesh );
}

/*! Surface meshing Delaunay refinement, only insert new vertices, no vertex is removed
*
*/
void _Surface_Delaunay_Mesh_Refinement( const char * _input_file, const char * _output_file )
{
    CDRTMesh dmesh;
	dmesh.read_m( _input_file );
	CDRTMesh tmesh;
	tmesh.Refine( dmesh );
	CDRTMesh omesh;
	tmesh.Convert( omesh );
	omesh.write_m( _output_file );
}


/*!
 *	Surface Delaunay remshing using Quad Tree data structure to improve the robustness
 */
void _Surface_Delaunay_Quad_Tree_Remesh( const char * _input_mesh, const char * _output_mesh )
{
		CQTDRMesh dmesh;
		dmesh.read_m( _input_mesh );


		//for( CQTDRMesh::MeshVertexIterator viter( &dmesh ); !viter.end(); viter ++ )
		//{
		//	CQTDelaunayRemeshVertex * pV = *viter;
		//	pV->uv() = pV->uv() * 2.0 - CPoint2(1,1);
		//}

		dmesh._lambda();
		dmesh._construct_tree( 8 );

		CQTDRMesh tmesh;
		tmesh.Remesh( dmesh );

		CQTDRMesh omesh;
		tmesh.Convert( omesh );
		omesh.write_m( _output_mesh );
}
