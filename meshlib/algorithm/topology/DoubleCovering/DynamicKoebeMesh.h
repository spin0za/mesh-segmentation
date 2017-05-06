#ifndef  _DYNAMIC_KOEBE_MESH_H_
#define  _DYNAMIC_KOEBE_MESH_H_

#include <map>
#include <vector>
#include <queue>

#include "Mesh/Vertex.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/boundary.h"
#include "Mesh/iterators.h"
#include "Mesh/BaseMesh.h"
#include "Mesh/DynamicMesh.h"


namespace MeshLib
{

namespace RicciFlow
{
/*! \brief CDynamicKoebeMesh class
*   
*  The mesh supports hole filling and hole puncture
*
*/
template<typename V, typename E, typename F, typename H>
class CDynamicKoebeMesh : public CDynamicMesh<V,E,F,H>
{
public:
	typedef V CVertex;
	typedef E CEdge;
	typedef F CFace;
	typedef H CHalfEdge;
	typedef CBoundary<V,E,F,H>  CBoundary;

	typedef MeshVertexIterator<V,E,F,H>        MeshVertexIterator;
	typedef MeshFaceIterator<V,E,F,H>          MeshFaceIterator;
	typedef MeshEdgeIterator<V,E,F,H>          MeshEdgeIterator;
	typedef FaceVertexIterator<V,E,F,H>        FaceVertexIterator;
	typedef FaceHalfedgeIterator<V,E,F,H>	   FaceHalfedgeIterator;
	typedef VertexEdgeIterator<V,E,F,H>		   VertexEdgeIterator;
	typedef VertexFaceIterator<V,E,F,H>		   VertexFaceIterator;
	typedef VertexVertexIterator<V,E,F,H>	   VertexVertexIterator;
	typedef VertexInHalfedgeIterator<V,E,F,H>	   VertexInHalfedgeIterator;
	typedef FaceEdgeIterator<V,E,F,H>		   FaceEdgeIterator;
	


	/*!
	 *	fill all the inner holes
	 */
	void _fill_holes();
	/*! 
	 * fill all the inner holes except the 0-th and the k-th
	 */
	void _fill_holes( int k );
	/*!
	 *	punch a hole at the first marker vertex v
	 */
	void _punch_holes( );
	/*!
	 *	punch all the holes 
	 */
	void _punch_all_holes( );
	/*!
	 *	read markers
	 */
	void _read_markers( const char * filename );
	/*! 
	 *	write markers
	 */
	void _rewrite_markers( const char * filename );
	void _update_markers( const char * filename );

	/*!
	 *	copy the position from the current mesh to the original mesh, 
	 *  \param pMesh - the original mesh with holes
	 */
	void _copy_position( CDynamicKoebeMesh<V,E,F,H> * pMesh );

protected:
	std::queue<V*> m_markers;
};

/*!
 *	embed the symmetric polygons along the symmetry axis
 */
template<typename V, typename E, typename F, typename H>
void CDynamicKoebeMesh<V,E,F,H>::_fill_holes()
{
	int vid = 0;
	for( MeshVertexIterator viter( this ); !viter.end(); viter ++ )
	{
		V * pv = *viter;
		vid = (vid > pv->id())?vid:pv->id();
	}

	int eid = 0;
	for( MeshEdgeIterator eiter( this ); !eiter.end(); eiter ++ )
	{
		E * pe = *eiter;
		eid = (eid > pe->id())?eid:pe->id();
	}

	int fid = 0;
	for( MeshFaceIterator fiter( this ); !fiter.end(); fiter ++ )
	{
		F * pf = *fiter;
		fid = (fid > pf->id())?fid:pf->id();
	}

	CBoundary boundary(this);

	  std::vector<CLoop<V,E,F,H>*> & loops = boundary.loops();

	  for( size_t i = 1; i < loops.size(); i ++ )
	  {
		  CLoop<V,E,F,H> * pL = loops[i];
			
		  std::vector<H*> hes;
		  CPoint sum(0,0,0);
		
		  for( std::list<H*>::iterator hiter = pL->halfedges().begin();
			  hiter != pL->halfedges().end(); hiter ++ )
		  {
			  H * pH = *hiter;
			  hes.push_back( pH );
			  V * pW = halfedgeTarget( pH );
			  sum += pW->point();
		  }
			
		  sum /= hes.size();

			// create vertex
			CVertex * v = new CVertex();
			assert( v != NULL );
			v->id() = ++vid;
			v->point() = sum;
			v->boundary() = false;
			v->halfedge() = NULL;

			m_verts.push_back( v );
			m_map_vert.insert( std::pair<int,CVertex*>(v->id(),v));
			m_markers.push( v );

			std::vector<E*> edges;
			for( size_t j = 0; j < hes.size(); j ++ )
			{

				// create edge 
				CEdge * e = new CEdge();	assert( e != NULL );
				e->id() = ++eid;
				m_edges.push_back( e );
				//m_map_edge.insert( std::pair<int,CEdge*>(e->id(),e));
				edges.push_back( e );

				// create halfedges & link to edge
				CHalfEdge * he[2];
				for (int k=0; k<2; k++)
				{
					he[k] = new CHalfEdge();	assert(he[k]);
					he[k]->vertex() = NULL;
					he[k]->edge() = e;
					he[k]->face() = NULL;
					e->halfedge(k) = he[k];
				}
			}

			std::vector<F*> faces;
			for( size_t j = 0; j < hes.size(); j ++ )
			{
				// create face & link he
				CFace * f = new CFace();
				assert( f != NULL );
				f->id() = ++fid;
				m_faces.push_back( f );
				m_map_face.insert( std::pair<int,tFace>(f->id(),f) );
				faces.push_back( f );
			}

			std::vector<H*> shs;
			for( size_t j = 0; j < hes.size(); j ++ )
			{
				// create face & link he
				CHalfEdge * ps = new CHalfEdge();
				assert( ps != NULL );
				CEdge     * pe = halfedgeEdge( hes[j] );
				pe->halfedge(1) = ps;
				ps->edge() = pe;
				ps->vertex() = hes[j]->he_prev()->vertex();
				shs.push_back( ps );
			}
			
			size_t n = hes.size();

			for( size_t j = 0; j < hes.size(); j ++ )
			{
				H * he[3];

				he[0] = edgeHalfedge(edges[j],1);
				he[1] = shs[j];
				he[2] = edgeHalfedge(edges[(j+1)%n],0);
			
				//halfedge->halfedge
				for( int k = 0; k < 3; k ++ )
				{
					he[k]->he_next() = he[(k+1)%3];
					he[k]->he_prev() = he[(k+2)%3];
				}
				//halfedge->face
				for( int k = 0; k < 3; k ++ )
				{
					he[k]->face() = faces[j];
				}
				faces[j]->halfedge() = he[0];

				//halfedge->vertex
				he[2]->vertex() = v;
				v->halfedge() = he[2];
				he[0]->vertex() = hes[j]->vertex();

			}

	  }	


}



/*!
 *	embed the symmetric polygons along the symmetry axis
 */
template<typename V, typename E, typename F, typename H>
void CDynamicKoebeMesh<V,E,F,H>::_punch_holes( )
{
	if( m_markers.empty()) 
	{
		std::cerr << "The markers are empty " << std::endl;
		return;
	}

	V * v = m_markers.front();
	m_markers.pop();

	std::vector<F*> faces;

	for( VertexFaceIterator vfiter( v ); !vfiter.end(); vfiter ++ )
	{
		F * pf = *vfiter;
		faces.push_back( pf );
	}

	std::vector<E*> edges;

	for( VertexEdgeIterator veiter( v ); !veiter.end(); veiter ++ )
	{
		E * pe = *veiter;
		edges.push_back( pe );
	}

	std::vector<H*> halfedges;

	for( size_t i = 0; i < faces.size(); i ++ )
	{
		F * pf = faces[i];

		for( FaceHalfedgeIterator fhiter( pf ); !fhiter.end(); fhiter ++)
		{
			H * ph = *fhiter;
			halfedges.push_back( ph );
		}
	}

	std::vector<H*> dual_hes;

	for( VertexInHalfedgeIterator vhiter( this,v ); !vhiter.end(); vhiter ++ )
	{
		H * ph = *vhiter;
		H * pw = (H*)ph->he_prev();
		dual_hes.push_back( pw );
	}


	for( size_t i = 0; i < faces.size(); i ++ )
	{
		F * pf = faces[i];
		m_faces.remove( pf );
		std::map<int,tFace>::iterator fiter = m_map_face.find( pf->id() );
		m_map_face.erase( fiter );
		delete pf;
	}

	for( size_t i = 0; i < edges.size(); i ++ )
	{
		E * pe = edges[i];
		m_edges.remove( pe );
		delete pe;
	}

	m_verts.remove( v );
	std::map<int,tVertex>::iterator viter = m_map_vert.find( v->id() );
	m_map_vert.erase( viter );
	delete v;

	for( size_t i = 0; i < dual_hes.size() ; i ++ )
	{
		H * ph = dual_hes[i];
		H * ps = halfedgeSym( ph );
		E * pe = halfedgeEdge( ph );
		pe->halfedge(0) = ps;
		pe->halfedge(1) = NULL;
	}

	for( size_t i = 0; i < halfedges.size(); i ++ )
	{
		H * ph = halfedges[i];
		delete ph;
	}
}


template<typename V, typename E, typename F, typename H>
void CDynamicKoebeMesh<V,E,F,H>::_read_markers( const char * input_file )
{
	std::fstream tis( input_file, std::fstream::in );

	if( tis.fail() )
	{
		std::cerr << "Error in opening file " << input_file << std::endl;
		return;
	}

	while( !tis.eof() )
	{
		int id;
		tis >> id;

		V * pV = idVertex( id );
		m_markers.push(pV);
	}

	tis.close();
};

template<typename V, typename E, typename F, typename H>
void CDynamicKoebeMesh<V,E,F,H>::_rewrite_markers( const char * output_file )
{
	std::fstream tos( output_file, std::fstream::out  );

	if( tos.fail() )
	{
		std::cerr << "Error in opening file " << output_file << std::endl;
		return;
	}

	while( !m_markers.empty())
	{
		V * pv = m_markers.front();
		m_markers.pop();
		tos << std::endl << pv->id();
		std::cout << pv->id() << std::endl;
	}

	tos.close();
};


template<typename V, typename E, typename F, typename H>
void CDynamicKoebeMesh<V,E,F,H>::_update_markers( const char * output_file )
{
	std::fstream tos( output_file, std::fstream::out  | std::fstream::app );

	if( tos.fail() )
	{
		std::cerr << "Error in opening file " << output_file << std::endl;
		return;
	}

	while( !m_markers.empty())
	{
		V * pv = m_markers.front();
		m_markers.pop();
		tos << std::endl << pv->id();
	}

	tos.close();
};

/*!
 *	embed the symmetric polygons along the symmetry axis
 */
template<typename V, typename E, typename F, typename H>
void CDynamicKoebeMesh<V,E,F,H>::_punch_all_holes( )
{
	while( !m_markers.empty() ) 
	{
		V * v = m_markers.front();
		m_markers.pop();

		std::vector<F*> faces;

		for( VertexFaceIterator vfiter( v ); !vfiter.end(); vfiter ++ )
		{
			F * pf = *vfiter;
			faces.push_back( pf );
		}

		std::vector<E*> edges;

		for( VertexEdgeIterator veiter( v ); !veiter.end(); veiter ++ )
		{
			E * pe = *veiter;
			edges.push_back( pe );
		}

		std::vector<H*> halfedges;

		for( size_t i = 0; i < faces.size(); i ++ )
		{
			F * pf = faces[i];

			for( FaceHalfedgeIterator fhiter( pf ); !fhiter.end(); fhiter ++)
			{
				H * ph = *fhiter;
				halfedges.push_back( ph );
			}
		}

		std::vector<H*> dual_hes;

		for( VertexInHalfedgeIterator vhiter( this,v ); !vhiter.end(); vhiter ++ )
		{
			H * ph = *vhiter;
			H * pw = (H*)ph->he_prev();
			dual_hes.push_back( pw );
		}


		for( size_t i = 0; i < faces.size(); i ++ )
		{
			F * pf = faces[i];
			m_faces.remove( pf );
			std::map<int,tFace>::iterator fiter = m_map_face.find( pf->id() );
			m_map_face.erase( fiter );
			delete pf;
		}

		for( size_t i = 0; i < edges.size(); i ++ )
		{
			E * pe = edges[i];
			m_edges.remove( pe );
			delete pe;
		}

		m_verts.remove( v );
		std::map<int,tVertex>::iterator viter = m_map_vert.find( v->id() );
		m_map_vert.erase( viter );
		delete v;

		for( size_t i = 0; i < dual_hes.size() ; i ++ )
		{
			H * ph = dual_hes[i];
			H * ps = halfedgeSym( ph );
			E * pe = halfedgeEdge( ph );
			pe->halfedge(0) = ps;
			pe->halfedge(1) = NULL;
		}

		for( size_t i = 0; i < halfedges.size(); i ++ )
		{
			H * ph = halfedges[i];
			delete ph;
		}
	}
}


/*!
 *	embed the symmetric polygons along the symmetry axis
 */
template<typename V, typename E, typename F, typename H>
void CDynamicKoebeMesh<V,E,F,H>::_copy_position( CDynamicKoebeMesh<V,E,F,H> * pMesh )
{
	for( MeshVertexIterator viter( pMesh ); !viter.end(); viter ++ )
	{
		V * pv = *viter;
		V * pw = idVertex( pv->id() );
		//std::cout << pv->id() << std::endl;
		pv->point() = pw->point();
	}
}


/*!
 *	embed the symmetric polygons along the symmetry axis
 */
template<typename V, typename E, typename F, typename H>
void CDynamicKoebeMesh<V,E,F,H>::_fill_holes( int k )
{
	int vid = 0;
	for( MeshVertexIterator viter( this ); !viter.end(); viter ++ )
	{
		V * pv = *viter;
		vid = (vid > pv->id())?vid:pv->id();
	}

	int eid = 0;
	for( MeshEdgeIterator eiter( this ); !eiter.end(); eiter ++ )
	{
		E * pe = *eiter;
		eid = (eid > pe->id())?eid:pe->id();
	}

	int fid = 0;
	for( MeshFaceIterator fiter( this ); !fiter.end(); fiter ++ )
	{
		F * pf = *fiter;
		fid = (fid > pf->id())?fid:pf->id();
	}

	 CBoundary boundary(this);

	  std::vector<CLoop<V,E,F,H>*> & loops = boundary.loops();

	  std::vector<CLoop<V,E,F,H>*> tloops;
	  tloops = loops;

	 for( size_t i = 0; i < tloops.size() - 1; i ++ )
	 {
		 for( size_t j = 0; j < tloops.size() - 1 - i; j ++ )
		 {
			 CLoop<V,E,F,H>* cL = tloops[j];
			 CLoop<V,E,F,H>* nL = tloops[j+1];
			
			 if( cL->halfedges().size() < nL->halfedges().size() )
			 {
				 tloops[j] = nL;
				 tloops[j+1] = cL;
			 }
		 }
	 }

	 for( size_t i = 0; i < tloops.size() - 1; i ++ )
	 {
		CLoop<V,E,F,H>* cL = tloops[i];
		std::cout << cL->halfedges().size() << std::endl;
	 }

	  for( size_t i = 1; i < tloops.size(); i ++ )
	  {
		  if( i == k ) continue;

		  CLoop<V,E,F,H> * pL = tloops[i];
			
		  std::vector<H*> hes;
		  CPoint sum(0,0,0);
		
		  for( std::list<H*>::iterator hiter = pL->halfedges().begin();
			  hiter != pL->halfedges().end(); hiter ++ )
		  {
			  H * pH = *hiter;
			  hes.push_back( pH );
			  V * pW = halfedgeTarget( pH );
			  sum += pW->point();
		  }
			
		  sum /= hes.size();

			// create vertex
			CVertex * v = new CVertex();
			assert( v != NULL );
			v->id() = ++vid;
			v->point() = sum;
			v->boundary() = false;
			v->halfedge() = NULL;

			m_verts.push_back( v );
			m_map_vert.insert( std::pair<int,CVertex*>(v->id(),v));
			m_markers.push( v );

			std::vector<E*> edges;
			for( size_t j = 0; j < hes.size(); j ++ )
			{

				// create edge 
				CEdge * e = new CEdge();	assert( e != NULL );
				e->id() = ++eid;
				m_edges.push_back( e );
				//m_map_edge.insert( std::pair<int,CEdge*>(e->id(),e));
				edges.push_back( e );

				// create halfedges & link to edge
				CHalfEdge * he[2];
				for (int k=0; k<2; k++)
				{
					he[k] = new CHalfEdge();	assert(he[k]);
					he[k]->vertex() = NULL;
					he[k]->edge() = e;
					he[k]->face() = NULL;
					e->halfedge(k) = he[k];
				}
			}

			std::vector<F*> faces;
			for( size_t j = 0; j < hes.size(); j ++ )
			{
				// create face & link he
				CFace * f = new CFace();
				assert( f != NULL );
				f->id() = ++fid;
				m_faces.push_back( f );
				m_map_face.insert( std::pair<int,tFace>(f->id(),f) );
				faces.push_back( f );
			}

			std::vector<H*> shs;
			for( size_t j = 0; j < hes.size(); j ++ )
			{
				// create face & link he
				CHalfEdge * ps = new CHalfEdge();
				assert( ps != NULL );
				CEdge     * pe = halfedgeEdge( hes[j] );
				pe->halfedge(1) = ps;
				ps->edge() = pe;
				ps->vertex() = hes[j]->he_prev()->vertex();
				shs.push_back( ps );
			}
			
			size_t n = hes.size();

			for( size_t j = 0; j < hes.size(); j ++ )
			{
				H * he[3];

				he[0] = edgeHalfedge(edges[j],1);
				he[1] = shs[j];
				he[2] = edgeHalfedge(edges[(j+1)%n],0);
			
				//halfedge->halfedge
				for( int k = 0; k < 3; k ++ )
				{
					he[k]->he_next() = he[(k+1)%3];
					he[k]->he_prev() = he[(k+2)%3];
				}
				//halfedge->face
				for( int k = 0; k < 3; k ++ )
				{
					he[k]->face() = faces[j];
				}
				faces[j]->halfedge() = he[0];

				//halfedge->vertex
				he[2]->vertex() = v;
				v->halfedge() = he[2];
				he[0]->vertex() = hes[j]->vertex();

			}

	  }	


}


typedef CDynamicKoebeMesh<CDynamicYamabeFlowVertex,CDynamicYamabeFlowEdge,CDynamicYamabeFlowFace,CDynamicYamabeFlowHalfEdge> CDKMesh;

} //namespace RicciFlow
} //namespace MeshLib

#endif