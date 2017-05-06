/*!
*      \file Stitch.cpp
*      \brief Implementation for CStitch class 
*	   \author David Gu
*      \date Documented on 10/12/2010
*
*	   Algorithm for stitching two meshes along their common boundaries.
*/

#include "Stitch.h"

using namespace MeshLib;

//CStitch constructor
//\param pMesh input base mesh
Topology::CStitch::CStitch( Topology::CSTMesh * pMesh ):m_pMesh( pMesh), m_boundary( m_pMesh )
{
}
//CStitch destructor
Topology::CStitch::~CStitch()
{
}

//planar area of the domain circled by a boundary loop hl

double Topology::CStitch::_area( std::list<Topology::CStitchHalfEdge*>  & hl )
{
	std::vector<CPoint2> points;

	double r = 0;
	int    n = 0;

	for( std::list<Topology::CStitchHalfEdge*>::iterator hiter = hl.begin(); hiter != hl.end(); hiter ++ )
	{
		Topology::CStitchHalfEdge * pH = *hiter;
		Topology::CStitchVertex   * pV = m_pMesh->halfedgeTarget( pH );

		CPoint2 p = pV->uv();
		points.push_back( p );
	}

	double s = 0;

	for( size_t i = 0; i < points.size(); i ++ )
	{
		CPoint2 B=points[i];
		CPoint2 C=points[(i+1)%points.size()];
		s += B[0] * C[1] - B[1] * C[0];
	}
	return s;
}

//convert a boundary loop,which is with greatest negative planar area, to a poly file

void Topology::CStitch::_hole_2_poly(  const char * output_name )
{

	std::vector<Topology::CSTMesh::CLoop*> & loops = m_boundary.loops();

	Topology::CSTMesh::CLoop * outer_loop = NULL;
	Topology::CSTMesh::CLoop * inner_loop = NULL;
	
	double min_area = 1e+10;

	for( size_t i = 0; i < loops.size(); i ++ )
	{
		Topology::CSTMesh::CLoop *pL = loops[i];
		double s = _area( pL->halfedges());
		printf("Area is %f\n", s );
		if( min_area > s )
		{ 
			min_area = s;
			inner_loop = pL; 
		};
	}

	std::vector<CPoint2> points;

	std::list<Topology::CStitchHalfEdge*>  & hl = inner_loop->halfedges();

	size_t n = 0;

	for( std::list<Topology::CStitchHalfEdge*>::iterator hiter = hl.begin(); hiter != hl.end(); hiter ++ )
	{
		Topology::CStitchHalfEdge * pH = *hiter;
		Topology::CStitchVertex   * pV = m_pMesh->halfedgeTarget( pH );
		points.push_back( pV->uv() );
		n ++;
	}

	std::fstream _os( output_name, std::fstream::out );
	if( _os.fail() )
	{
		fprintf(stderr,"Error is opening file %s\n", output_name );
		return;
	}
	_os << points.size() << " 2 0 0" << std::endl;

	//FILE * fp = fopen( output_name, "w");
	//fprintf(fp,"%d 2 0 0\n", points.size());
	
	for( size_t i = 0; i < points.size(); i ++ )
	{
		//fprintf(fp, "%d %g %g\n", i+1, points[i][0], points[i][1] );
		_os << i+1 << " " << points[i][0] << " " << points[i][1] << std::endl;
	}

	//fprintf(fp, "%d 0\n", n);
	
	_os << n << " 0" << std::endl;

	for( size_t i = 0; i < n-1; i ++ )
	{
		//fprintf(fp, "%d %d %d\n", i+1, i+1, i + 2 );

		_os << i+1 << " " << i + 1 << " " << i + 2 << std::endl;
	}
	
	_os << n << " " << n << " 1" << std::endl;

	//fprintf(fp, "%d %d %d\n", n, n, 1 );
	//fprintf(fp, "0\n" );
	
	_os << "0" << std::endl;
	_os.close();
	
	//fclose(fp);

}

//stitch a patch mesh patch to the base mesh, output to the merged mesh
// patch is the patch mesh
// output_name is the name of the output mesh

void Topology::CStitch::_stitch( Topology::CSTMesh & patch, const char * output_name )
{

	std::vector<Topology::CSTMesh::CLoop*> & loops = m_boundary.loops();

	Topology::CSTMesh::CLoop * outer_loop = NULL;
	Topology::CSTMesh::CLoop * inner_loop = NULL;
	
	double min_area = 1e+10;

	for( size_t i = 0; i < loops.size(); i++ )
	{
		Topology::CSTMesh::CLoop *pL = loops[i];
		double s = _area( pL->halfedges());
		printf("Area is %f\n", s );
		if( min_area > s )
		{ 
			min_area = s;
			inner_loop = pL; 
		};
	}

	std::vector<Topology::CStitchVertex*> common_verts;


	std::list<Topology::CStitchHalfEdge*>  & hl = inner_loop->halfedges();

	for( std::list<Topology::CStitchHalfEdge*>::iterator hiter = hl.begin(); hiter != hl.end(); hiter ++ )
	{
		Topology::CStitchHalfEdge * pH = *hiter;
		Topology::CStitchVertex   * pV = m_pMesh->halfedgeTarget( pH );
		common_verts.push_back( pV );
	}

	int max_vid = -1;

	for( Topology::CSTMesh::MeshVertexIterator	viter( m_pMesh ); !viter.end(); viter ++ )
	{
		Topology::CStitchVertex * pV = *viter;
		int id = pV->id();
		pV->idx() = id;

		max_vid = (max_vid > id )?max_vid:id;
	}

	int n = common_verts.size();

	for( Topology::CSTMesh::MeshVertexIterator	viter( &patch); !viter.end(); viter ++ )
	{
		Topology::CStitchVertex * pV = *viter;
	
		int id = pV->id();
		if( id <= n )
		{
			Topology::CStitchVertex * pMatch = common_verts[id-1];
			pV->idx() = pMatch->idx();
		}
		else
		{
			pV->idx()  = id - n + max_vid;
		}
	}


	CSTMesh  union_mesh;


	for( CSTMesh::MeshVertexIterator	viter( m_pMesh ); !viter.end(); viter ++ )
	{
		Topology::CStitchVertex * pV = *viter;
		int id = pV->idx();
		
		Topology::CStitchVertex * pW = union_mesh.createVertex( id );
		pW->point()  = pV->point();
		pW->string() = pV->string();
	}


	for( CSTMesh::MeshVertexIterator	viter( &patch); !viter.end(); viter ++ )
	{
		Topology::CStitchVertex * pV = *viter;
		int id = pV->id();
		CPoint p = pV->point();
		if( id <= n ) continue;
		id = pV->idx();

		Topology::CStitchVertex * pW = union_mesh.createVertex( id );
		pW->point()  = pV->point();
		pW->string() = pV->string();

	}


	int max_fid = -1;

	for( CSTMesh::MeshFaceIterator	fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		CStitchFace * f = *fiter;
		int fid = f->id();

		Topology::CStitchVertex * w[3];
		int i = 0;
		for( CSTMesh::FaceVertexIterator fviter( f ); !fviter.end(); fviter ++ )
		{
			Topology::CStitchVertex * pV = *fviter;
			int id = pV->idx();
			Topology::CStitchVertex * pW = union_mesh.idVertex( id );
			w[i++] = pW;
		}
		
		f->idx() = fid;
		CStitchFace * pF = union_mesh.createFace( w, fid );
	
		max_fid = (max_fid > fid )?max_fid:fid;
	}

	for( CSTMesh::MeshFaceIterator	fiter( &patch); !fiter.end(); fiter ++ )
	{
		CStitchFace * f = *fiter;
		int fid = f->id() + max_fid;

		Topology::CStitchVertex * w[3];
		int i = 0;
		for( CSTMesh::FaceVertexIterator fviter( f ); !fviter.end(); fviter ++ )
		{
			Topology::CStitchVertex * pV = *fviter;
			int id = pV->idx();
			Topology::CStitchVertex * pW = union_mesh.idVertex( id );
			w[i++] = pW;
		}
		
		f->idx() = fid;
		CStitchFace * pF = union_mesh.createFace( w, fid );
	}

	//copy edge string
	for( CSTMesh::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		Topology::CStitchEdge * e = *eiter;

		Topology::CStitchVertex   * v1 = m_pMesh->edgeVertex1( e );
		Topology::CStitchVertex   * v2 = m_pMesh->edgeVertex2( e );

		int idx1 = v1->idx();
		int idx2 = v2->idx();

		Topology::CStitchVertex * w1 = union_mesh.idVertex( idx1 );
		Topology::CStitchVertex * w2 = union_mesh.idVertex( idx2 );

		Topology::CStitchEdge * pWE = union_mesh.vertexEdge( w1, w2 );

		pWE->string() = e->string();
	}


	for( CSTMesh::MeshEdgeIterator eiter( &patch ); !eiter.end(); eiter ++ )
	{
		Topology::CStitchEdge * e = *eiter;

		Topology::CStitchVertex   * v1 = m_pMesh->edgeVertex1( e );
		Topology::CStitchVertex   * v2 = m_pMesh->edgeVertex2( e );

		int idx1 = v1->idx();
		int idx2 = v2->idx();

		Topology::CStitchVertex * w1 = union_mesh.idVertex( idx1 );
		Topology::CStitchVertex * w2 = union_mesh.idVertex( idx2 );

		Topology::CStitchEdge * pWE = union_mesh.vertexEdge( w1, w2 );

		pWE->string() = e->string();
	}


	for( CSTMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		Topology::CStitchVertex *pV  = *viter;
		int vid = pV->idx();
		Topology::CStitchVertex * pW = union_mesh.idVertex( vid ); 
		pW->string() = pV->string();
	}


	for( CSTMesh::MeshVertexIterator viter( &patch ); !viter.end(); viter ++ )
	{
		Topology::CStitchVertex *pV  = *viter;
		int vid = pV->idx();
		Topology::CStitchVertex * pW = union_mesh.idVertex( vid ); 
		pW->string() = pV->string();
	}


	for( std::list<Topology::CStitchHalfEdge*>::iterator hiter = hl.begin(); hiter != hl.end(); hiter ++ )
	{
		Topology::CStitchHalfEdge * pH = *hiter;
		Topology::CStitchEdge     * edge = m_pMesh->halfedgeEdge( pH );
		
		Topology::CStitchVertex   * v1 = m_pMesh->edgeVertex1( edge );
		Topology::CStitchVertex   * v2 = m_pMesh->edgeVertex2( edge );
		
		int idx1 = v1->idx();
		int idx2 = v2->idx();
		
		Topology::CStitchVertex  * w1 = union_mesh.idVertex( idx1 );
		Topology::CStitchVertex  * w2 = union_mesh.idVertex( idx2 );

		Topology::CStitchEdge * pWE = union_mesh.vertexEdge( w1, w2 );
		pWE->sharp() = true;
	}


	//copy corner angle

	for( CSTMesh::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		CStitchFace * pF = *fiter;
		CStitchFace * pWF = union_mesh.idFace( pF->idx() );

		for( CSTMesh::FaceVertexIterator fviter( pF ); !fviter.end(); fviter ++ )
		{
			Topology::CStitchVertex * pV = * fviter;
			Topology::CStitchVertex * pW = union_mesh.idVertex( pV->idx() );
			Topology::CStitchHalfEdge *  pH = m_pMesh->corner( pV, pF );
			Topology::CStitchHalfEdge * pWH = union_mesh.corner( pW, pWF );
			pWH->string() = pH->string();
			pWH->angle() = pH->angle();
		}
	}

	for( CSTMesh::MeshFaceIterator fiter( &patch ); !fiter.end(); fiter ++ )
	{
		CStitchFace * pF = *fiter;
		CStitchFace * pWF = union_mesh.idFace( pF->idx() );

		for( CSTMesh::FaceVertexIterator fviter( pF ); !fviter.end(); fviter ++ )
		{
			Topology::CStitchVertex * pV = * fviter;
			Topology::CStitchVertex * pW = union_mesh.idVertex( pV->idx() );
			Topology::CStitchHalfEdge *  pH = patch.corner( pV, pF );
			Topology::CStitchHalfEdge * pWH = union_mesh.corner( pW, pWF );

			pWH->string() = pH->string();
			pWH->angle() = pH->angle();
		}
	}

	union_mesh.write_m( output_name );
}


//convert a boundary loop,which is with greatest negative planar area, to a poly file

void Topology::CStitch::_fill_the_biggest_hole(  const char * output_name )
{

	std::vector<Topology::CSTMesh::CLoop*> & loops = m_boundary.loops();

	Topology::CSTMesh::CLoop * outer_loop = NULL;
	Topology::CSTMesh::CLoop * inner_loop = NULL;
	
	double min_area = 1e+10;

	for( size_t i = 0; i < loops.size(); i ++ )
	{
		Topology::CSTMesh::CLoop *pL = loops[i];
		double s = _area( pL->halfedges());
		printf("Area is %f\n", s );
		if( min_area > s )
		{ 
			min_area = s;
			inner_loop = pL; 
		};
	}

	std::vector<CPoint2> points;
	std::vector<Topology::CStitchVertex*> vertices;

	std::list<Topology::CStitchHalfEdge*>  & hl = inner_loop->halfedges();

	size_t n = 0;

	for( std::list<Topology::CStitchHalfEdge*>::iterator hiter = hl.begin(); hiter != hl.end(); hiter ++ )
	{
		Topology::CStitchHalfEdge * pH = *hiter;
		Topology::CStitchVertex   * pV = m_pMesh->halfedgeTarget( pH );
		vertices.push_back( pV );
		points.push_back( pV->uv() );
		n ++;
	}
	
	CPoly poly;

	for( size_t i = 0; i < points.size(); i ++ )
	{
		poly.points().push_back( points[i] );
	}

	for( size_t i = 0; i < n-1; i ++ )
	{
		poly.segments().push_back( CPolySegment(i+1,i+2) );
	}
	poly.segments().push_back( CPolySegment(n,1) );
	
	CDTMesh Delaunay;
	CDTMesh patch;

	Delaunay.ProcessPoly( poly, patch );
	//patch.write_m("debug.m");
	
	for( CDTMesh::MeshVertexIterator viter( &patch ); !viter.end(); viter ++ )
	{
		CDelaunayVertex * pV = *viter;
		if( pV->boundary() )
		{
			pV->father() = vertices[pV->id()-4]->id();
		}
		else
		{
			pV->father() = -1;
		}
	}

	

	CSTMesh union_mesh;
	_union( &patch, union_mesh );

	union_mesh.write_m( output_name );
}


//stitch a patch mesh patch to the base mesh, output to the merged mesh
// patch is the patch mesh
// output_name is the name of the output mesh

void Topology::CStitch::_union( CDTMesh *  patch, Topology::CSTMesh & union_mesh )
{

	int n  = 0;
	for( CDTMesh::MeshVertexIterator viter( patch ); !viter.end(); viter ++ )
	{
		CDelaunayVertex * pV = *viter;
		if( pV->boundary() ) n++;
	}

	int max_vid = -1;
	for( Topology::CSTMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		Topology::CStitchVertex * pV = *viter;
		max_vid = (max_vid > pV->id() )?max_vid:pV->id();
	}

	for( CDTMesh::MeshVertexIterator viter( patch); !viter.end(); viter ++ )
	{
		CDelaunayVertex * pV = *viter;
		if( pV->boundary() )
		{
			printf("%d - %d\n", pV->id(), pV->father() );
			pV->id() = pV->father();
		}
		else
		{
			pV->id()  = pV->id() - n + max_vid;
		}
	}

	for( Topology::CSTMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		Topology::CStitchVertex * pV = *viter;
		int id = pV->id();
		
		Topology::CStitchVertex * pW = union_mesh.createVertex( id );
		//pW->point()  = pV->point();
		pW->point()  = CPoint(pV->uv()[0], pV->uv()[1],0);
		pW->rgb() = pV->rgb();
		pW->string() = pV->string();
	}


	for( CDTMesh::MeshVertexIterator	viter( patch); !viter.end(); viter ++ )
	{
		CDelaunayVertex * pV = *viter;
		if( pV->boundary() ) continue;
		int id = pV->id();
		CPoint p = pV->point();

		Topology::CStitchVertex * pW = union_mesh.createVertex( id );
		//pW->point()  = pV->point();
		pW->point()  = CPoint(pV->uv()[0], pV->uv()[1],0);
		pW->string() = pV->string();
	}


	int max_fid = -1;
	int max_sid = -1;

	for( Topology::CSTMesh::MeshFaceIterator	fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		Topology::CStitchFace * f = *fiter;
		int fid = f->id();
		int sid = f->segment();
		max_sid = ( max_sid > sid )? max_sid:sid;	
		max_fid = ( max_fid > fid )? max_fid:fid;
	}

	for( Topology::CSTMesh::MeshFaceIterator	fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		Topology::CStitchFace * f = *fiter;
		int fid = f->id();

		Topology::CStitchVertex * w[3];
		int i = 0;
		for( Topology::CSTMesh::FaceVertexIterator fviter( f ); !fviter.end(); fviter ++ )
		{
			Topology::CStitchVertex * pV = *fviter;
			int id = pV->id();
			Topology::CStitchVertex * pW = union_mesh.idVertex( id );
			w[i++] = pW;
		}
		
		Topology::CStitchFace * pF = union_mesh.createFace( w, fid );
		pF->segment() = f->segment();
		pF->string()  = f->string();
	}
	double rgbs[6][3]={{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1},{0,1,1}};

	for( CDTMesh::MeshFaceIterator	fiter( patch); !fiter.end(); fiter ++ )
	{
		CDelaunayFace * f = *fiter;
		int fid = f->id() + max_fid;
		f->id() = fid;

		Topology::CStitchVertex * w[3];
		int i = 0;
		for( CDTMesh::FaceVertexIterator fviter( f ); !fviter.end(); fviter ++ )
		{
			CDelaunayVertex * pV = *fviter;
			int id = pV->id();
			Topology::CStitchVertex * pW = union_mesh.idVertex( id );
			w[i++] = pW;
			if( pV->boundary() ) continue;
			int c = (max_sid+1)%6;
			pW->rgb() = CPoint(rgbs[c][0], rgbs[c][1], rgbs[c][2]);
		}
		
		Topology::CStitchFace * pF = union_mesh.createFace( w, fid );
		pF->segment() = max_sid + 1;
		pF->string() = f->string();
	}

	//copy edge string
	for( Topology::CSTMesh::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		Topology::CStitchEdge * e = *eiter;

		Topology::CStitchVertex   * v1 = m_pMesh->edgeVertex1( e );
		Topology::CStitchVertex   * v2 = m_pMesh->edgeVertex2( e );

		int idx1 = v1->id();
		int idx2 = v2->id();

		Topology::CStitchVertex * w1 = union_mesh.idVertex( idx1 );
		Topology::CStitchVertex * w2 = union_mesh.idVertex( idx2 );

		Topology::CStitchEdge * pWE = union_mesh.vertexEdge( w1, w2 );

		pWE->string() = e->string();
	}


	for( CDTMesh::MeshEdgeIterator eiter( patch ); !eiter.end(); eiter ++ )
	{
		CDelaunayEdge * e = *eiter;

		CDelaunayVertex   * v1 = patch->edgeVertex1( e );
		CDelaunayVertex   * v2 = patch->edgeVertex2( e );

		int idx1 = v1->id();
		int idx2 = v2->id();

		Topology::CStitchVertex * w1 = union_mesh.idVertex( idx1 );
		Topology::CStitchVertex * w2 = union_mesh.idVertex( idx2 );

		Topology::CStitchEdge * pWE = union_mesh.vertexEdge( w1, w2 );

		pWE->string() = e->string();
	}


	for( Topology::CSTMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		Topology::CStitchVertex *pV  = *viter;
		int vid = pV->id();
		Topology::CStitchVertex * pW = union_mesh.idVertex( vid ); 
		pW->string() = pV->string();
	}


	for( CDTMesh::MeshVertexIterator viter( patch ); !viter.end(); viter ++ )
	{
		CDelaunayVertex *pV  = *viter;
		int vid = pV->id();
		Topology::CStitchVertex * pW = union_mesh.idVertex( vid ); 
		pW->string() = pV->string();
	}

	for( CDTMesh::MeshEdgeIterator eiter( patch ); !eiter.end(); eiter ++ )
	{

		CDelaunayEdge * pE = *eiter;
		if( !pE->boundary() ) continue;

		CDelaunayVertex   * v1 = patch->edgeVertex1( pE );
		CDelaunayVertex   * v2 = patch->edgeVertex2( pE );
		
		int id1 = v1->id();
		int id2 = v2->id();
		
		Topology::CStitchVertex  * w1 = union_mesh.idVertex( id1 );
		Topology::CStitchVertex  * w2 = union_mesh.idVertex( id2 );

		Topology::CStitchEdge * pWE = union_mesh.vertexEdge( w1, w2 );
		pWE->sharp() = true;
	}


	//copy corner angle

	for( Topology::CSTMesh::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		Topology::CStitchFace * pF = *fiter;
		Topology::CStitchFace * pWF = union_mesh.idFace( pF->id() );

		for( Topology::CSTMesh::FaceVertexIterator fviter( pF ); !fviter.end(); fviter ++ )
		{
			Topology::CStitchVertex * pV = * fviter;
			Topology::CStitchVertex * pW = union_mesh.idVertex( pV->id() );
			Topology::CStitchHalfEdge *  pH = m_pMesh->corner( pV, pF );
			Topology::CStitchHalfEdge * pWH = union_mesh.corner( pW, pWF );
			pWH->string() = pH->string();
			pWH->angle() = pH->angle();
		}
	}

	for( CDTMesh::MeshFaceIterator fiter( patch ); !fiter.end(); fiter ++ )
	{
		CDelaunayFace * pF = *fiter;
		Topology::CStitchFace * pWF = union_mesh.idFace( pF->id() );

		for( CDTMesh::FaceVertexIterator fviter( pF ); !fviter.end(); fviter ++ )
		{
			CDelaunayVertex * pV = * fviter;
			Topology::CStitchVertex * pW = union_mesh.idVertex( pV->id() );
			CDelaunayHalfEdge *  pH = patch->corner( pV, pF );
			Topology::CStitchHalfEdge * pWH = union_mesh.corner( pW, pWF );

			pWH->string() = pH->string();
			//pWH->angle() = pH->angle();
		}
	}
}

//punch hole, remove all faces with segment id
void Topology::CStitch::_remove_segment( int segment_id, const char * output_name )
{

	for( Topology::CSTMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		Topology::CStitchVertex * pV = *viter;
		pV->touched() = false;
	}

	for( Topology::CSTMesh::MeshFaceIterator	fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		Topology::CStitchFace * f = *fiter;
		if( f->segment() == segment_id ) continue;

		for( Topology::CSTMesh::FaceVertexIterator fviter( f ); !fviter.end(); fviter ++ )
		{
			Topology::CStitchVertex * pV = *fviter;
			pV->touched() = true;
		}
	}

	Topology::CSTMesh mesh;

	for( Topology::CSTMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		Topology::CStitchVertex * pV = *viter;
		if( !pV->touched() ) continue;
		int id = pV->id();
		
		Topology::CStitchVertex * pW = mesh.createVertex( id );
		pW->point()  = pV->point();
		pW->string() = pV->string();
		pW->rgb()    = pV->rgb();
		pW->uv()     = pV->uv();
	}


	for( Topology::CSTMesh::MeshFaceIterator	fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		Topology::CStitchFace * f = *fiter;
		if( f->segment() == segment_id  ) continue;
		int fid = f->id();

		Topology::CStitchVertex * w[3];
		int i = 0;
		for( Topology::CSTMesh::FaceVertexIterator fviter( f ); !fviter.end(); fviter ++ )
		{
			Topology::CStitchVertex * pV = *fviter;
			int id = pV->id();
			Topology::CStitchVertex * pW = mesh.idVertex( id );
			w[i++] = pW;
		}
		
		Topology::CStitchFace * pF = mesh.createFace( w, fid );
		pF->segment() = f->segment();
	}

	//copy edge string
	for( Topology::CSTMesh::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		Topology::CStitchEdge * e = *eiter;

		Topology::CStitchFace * f1 = m_pMesh->edgeFace1( e );
		Topology::CStitchFace * f2 = m_pMesh->edgeFace1( e );
		if( f2 == NULL && f1->segment() == segment_id ) continue;
		if( f2 != NULL && f1->segment() == segment_id && f2->segment() == segment_id ) continue;

		Topology::CStitchVertex   * v1 = m_pMesh->edgeVertex1( e );
		Topology::CStitchVertex   * v2 = m_pMesh->edgeVertex2( e );

		int idx1 = v1->id();
		int idx2 = v2->id();

		Topology::CStitchVertex * w1 = mesh.idVertex( idx1 );
		Topology::CStitchVertex * w2 = mesh.idVertex( idx2 );
		Topology::CStitchEdge   * pWE = mesh.vertexEdge( w1, w2 );

		pWE->string() = e->string();
	}






	//copy corner angle

	for( Topology::CSTMesh::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		Topology::CStitchFace * pF = *fiter;
		if( pF->segment() == segment_id ) continue;

		Topology::CStitchFace * pWF = mesh.idFace( pF->id() );

		for( Topology::CSTMesh::FaceVertexIterator fviter( pF ); !fviter.end(); fviter ++ )
		{
			Topology::CStitchVertex * pV = * fviter;
			Topology::CStitchVertex * pW = mesh.idVertex( pV->id() );
			Topology::CStitchHalfEdge *  pH = m_pMesh->corner( pV, pF );
			Topology::CStitchHalfEdge * pWH = mesh.corner( pW, pWF );
			pWH->string() = pH->string();
			pWH->angle() = pH->angle();
		}
	}

	mesh.write_m( output_name );
}

//Copy UV
void Topology::CStitch::_copy_uv( CSTMesh * mesh )
{
	for( Topology::CSTMesh::MeshVertexIterator viter( mesh ); !viter.end(); viter ++ )
	{
		Topology::CStitchVertex * pV = *viter;
		Topology::CStitchVertex * pW = m_pMesh->idVertex( pV->id() );
		CPoint pos = pW->point();
		pV->uv() = CPoint2( pos[0], pos[1] );

		CParser parser( pV->string() );
		parser._removeToken( "uv" );
		parser._toString( pV->string() );

		std::stringstream iss;
		iss << "uv=(" << pV->uv()[0]<< " "<< pV->uv()[1] << ")";

		if( pV->string().length() > 0 )
		{
			pV->string() += " ";
		}
		pV->string() += iss.str();
	}
}

