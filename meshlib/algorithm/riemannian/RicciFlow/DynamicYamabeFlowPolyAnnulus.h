/*! \file   DynamicYamabeFlowPolyAnnulus.h
 *  \brief  Dynamic Euclidean Yamabe flow for poly annulus
 *  \author David Gu
 *  \date   documented on 10/20/2013
 *
 *	Algorithm for Dynamic Yamabe Flow for Poly Annulus
 */

#ifndef _DYNAMIC_YAMABE_FLOW_POLY_ANNULUS_H_
#define _DYNAMIC_YAMABE_FLOW_POLY_ANNULUS_H_

#include "Topology/DoubleCovering/DynamicDoubleCoveringMesh.h"
#include "Topology/DoubleCovering/DynamicDivideMesh.h"
#include "DynamicYamabeFlow.h"


namespace MeshLib
{
namespace RicciFlow
{

template<typename V>
class CEntry
{
public:
	CEntry( V* pointer, int father ){ m_point = pointer; m_father = father; };
	~CEntry(){};
public:
	V *    m_point;
	int    m_father;
	bool operator<(const CEntry & e ) const
	{
		return m_father < e.m_father;
	};
	bool operator>(const CEntry & e ) const 
	{
		return m_father > e.m_father;
	};

};


	/*! \brief CDynamicYamabeFlowExtremalLength class
	 *  
	 *  Algorithm for Dynamic Euclidean Yamabe flow for Extremal Length
	 */
	 
 template<typename M>
 class CDynamicYamabeFlowPolyAnnulus : public CDynamicYamabeFlow<M>
  {
  public:
	  typename M::CVertex *  TV;

	  /*! \brief CDynamicYamabeFlowPolyAnnulus constructor
	   *
	   *  call base class constructor 
	   */
	  CDynamicYamabeFlowPolyAnnulus( M * p_half_Mesh, M * p_double_Mesh );
	
	  /*! compute the metric to circle domain on a cylinder 
	  */
	  void _calculate_metric_circle_domain_on_cylinder();

  protected:


	/*!
	 *	Set the target curvature on each vertex
	 */
    virtual void    _set_target_curvature();

	M * m_pHalfMesh;
	typename M::CBoundary m_boundary;

  private:
	  //two vertices are dual to each other, if they share the same father
	  //the original boundary vertices on the half mesh are self-dual.
	  void _match_dual_vertex();
	  //each boundry vertex on half mesh has an ancestor on the double mesh
	  //each self-dual vertex on the double mesh has an ancestor on the half mesh
	  void _match_ancestor_vertex();
	  //match the ancestor edge
	  typename M::CEdge * _match_ancestor_edge( typename M::CHalfEdge * pH );

	  //set the target curvature on the half mesh	
	  void _set_half_mesh_target_curvature();
	  
	  //copy target curvature from the half mesh to the double mesh
	  void _copy_target_curvature_from_half_to_double();
	  //copy edge length from the double mesh to the half mesh
	  void _copy_edge_length_from_double_to_half();

  private:
	  //edges cross the boundary
	  std::set<typename M::CEdge*> m_cross_edges;
	  void _locate_cross_edges();
	  typename M::CVertex *  _corner_edge(typename M::CEdge * e );
	  void _trace_swapped_edge( typename M::CHalfEdge * pH );
	  
	  double _trace( typename M::CEdge * e );

 private:
	 typename M::CHalfEdge * _parallel_transport_upper( typename M::CHalfEdge * h )
	{
		return m_pMesh->halfedgeNext( m_pMesh->halfedgeSym( m_pMesh->halfedgePrev( m_pMesh->halfedgeSym(h) ) ) );
	};
	 typename M::CHalfEdge * _parallel_transport_lower( typename M::CHalfEdge * h )
	{
		return m_pMesh->halfedgePrev( m_pMesh->halfedgeSym( m_pMesh->halfedgeNext( m_pMesh->halfedgeSym(h) ) ) );
	};
	 typename M::CVertex * _third_vertex( typename M::CHalfEdge * h )
	{
		return m_pMesh->halfedgeTarget( m_pMesh->halfedgeNext(h) );
	};

private:
	//embedding a sequence of ajdacent faces

	//embed the first face
	void _embed_first_face( typename M::CFace * head );
	//embed other face
	void _embed_face( typename M::CFace * head );	
	/*!
	 *	embed a symmetric polygon
	 */
	void _embed_symmetric_polygon( std::vector<typename M::CHalfEdge *> & half_edges );

 };


/*!
 *	an edge is a corner edge if its 3rd vertex is self-dual
 */
template<typename M>
typename M::CVertex * CDynamicYamabeFlowPolyAnnulus<M>::_corner_edge( typename M::CEdge * e )
{
	M::CHalfEdge * h = m_pMesh->edgeHalfedge( e, 0 );
	M::CVertex   * v = _third_vertex( h );
	if( v->dual() == v ) return v;

	h = m_pMesh->edgeHalfedge( e, 1 );
	v = _third_vertex( h );
	if( v->dual() == v ) return v;
	
	return NULL;
};

template<typename M>
void CDynamicYamabeFlowPolyAnnulus<M>::_locate_cross_edges()
{
	m_cross_edges.clear();

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge   * e = *eiter;
		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );
		if( v1 == v2 ) continue;

		if( v1->dual() == v2 && v2->dual() == v1 )
			m_cross_edges.insert( e );
	}
};

template<typename M>
CDynamicYamabeFlowPolyAnnulus<M>::CDynamicYamabeFlowPolyAnnulus( M * p_half_Mesh, M * p_double_Mesh ):m_boundary( p_half_Mesh ), CDynamicYamabeFlow<M>( p_double_Mesh )
{ 
	m_pHalfMesh = p_half_Mesh; 
	_match_dual_vertex();
	_match_ancestor_vertex();
	
	//initialize edge length

	for( M::MeshEdgeIterator eiter( m_pHalfMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		e->length() = m_pHalfMesh->edgeLength( e );
	}

	//initialize edge length

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		e->length() = m_pHalfMesh->edgeLength( e );
	}

};


//set target curvature

template<typename M>
void CDynamicYamabeFlowPolyAnnulus<M>::_set_half_mesh_target_curvature()
{
  for( M::MeshVertexIterator viter( m_pHalfMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
	  v->target_k() = 0;
  }


  std::vector<M::CLoop*> & loops = m_boundary.loops();

  int i = 0;

  for( size_t i = 0; i < loops.size(); i ++ )
  {
	  if( i < 2 ) continue;
	  M::CLoop * pL = loops[i];
		
	  std::vector<M::CHalfEdge*> hes;

	  for( std::list<M::CHalfEdge*>::iterator hiter = pL->halfedges().begin();
		  hiter != pL->halfedges().end(); hiter ++ )
	  {
		  M::CHalfEdge * pH = *hiter;
		  hes.push_back( pH );
	  }
		
	  double length = 0;
	  for( size_t i = 0; i < hes.size(); i ++ )
	  {
		  M::CHalfEdge * pH = hes[i];
		  M::CEdge * e = m_pHalfMesh->halfedgeEdge( pH );
		  std::cout << e->length() << std::endl;
		  length += e->length();
	  }
	
	  for( size_t i = 0; i < hes.size(); i ++ )
	  {
		  M::CHalfEdge * pH = hes[i];
		  M::CHalfEdge * nH = hes[(i+1)%hes.size()];
		  M::CVertex   * pV = m_pHalfMesh->halfedgeTarget( pH );

		  M::CEdge * pE = m_pHalfMesh->halfedgeEdge( pH );
		  M::CEdge * nE = m_pHalfMesh->halfedgeEdge( nH );

		  pV->target_k() = -2 * PI * ( pE->length() + nE->length() )/( 2 * length );
	  }


  }
};


//set target curvature

template<typename M>
void CDynamicYamabeFlowPolyAnnulus<M>::_copy_target_curvature_from_half_to_double()
{

  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
	  M::CVertex * w = m_pHalfMesh->idVertex( v->father() );
	if( m_pHalfMesh->isBoundary( w ) )
	{
		v->target_k() = w->target_k() * 2.0;
	}
	else
	{
		v->target_k() = w->target_k();
	}
  }

};

//set target curvature

template<typename M>
void CDynamicYamabeFlowPolyAnnulus<M>::_set_target_curvature()
{

  _copy_edge_length_from_double_to_half();
  _set_half_mesh_target_curvature();
  _copy_target_curvature_from_half_to_double();
};

/*!
 *	two vertices are dual to each other, if they share the same father,
 *	Assume the father field has been set
 *  The vertex dual field will be set
 */

template<typename M>
void CDynamicYamabeFlowPolyAnnulus<M>::_match_dual_vertex()
{
	std::priority_queue<CEntry<M::CVertex> > pq;

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		std::cout << pV->id() << std::endl;
		CEntry<M::CVertex> entry(pV, pV->father());
		pq.push( entry );
	}
	std::vector<M::CVertex*> vs;
	while (! pq.empty()) {
		CEntry<M::CVertex> e = pq.top();
		M::CVertex * pV =  e.m_point;
		vs.push_back( pV );
		pq.pop();
	}
	for( size_t i = 0; i < vs.size(); i ++ )
	{
		M::CVertex * pV = vs[i];
		M::CVertex * pW = vs[(i+1)%vs.size()];

		if( pV->father() != pW->father() )
		{
			pV->dual() = pV;
		}
		else
		{
			pV->dual() = pW;
			pW->dual() = pV;
			i ++;
		}
	}

};

//set target curvature

template<typename M>
void CDynamicYamabeFlowPolyAnnulus<M>::_match_ancestor_vertex()
{
	//uniformization metric
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
	  M::CVertex * w = m_pHalfMesh->idVertex( v->father() );
	if( m_pHalfMesh->isBoundary( w ) )
	{
		v->ancestor() = w;
		w->ancestor() = v;
	}
  }

};


/*!
 *	two vertices are dual to each other, if they share the same father,
 *	Assume the father field has been set
 *  The vertex dual field will be set
 */

template<typename M>
void CDynamicYamabeFlowPolyAnnulus<M>::_copy_edge_length_from_double_to_half()
{
  //find all cross edges	
  _locate_cross_edges();

  std::vector<M::CLoop*> & loops = m_boundary.loops();

  int i = 0;

  for( size_t i = 0; i < loops.size(); i ++ )
  {
	  if( i < 2 ) continue;
	  M::CLoop * pL = loops[i];
	
	  //contains all swapped edges
	  std::vector<M::CHalfEdge*> swapped_hes;

	  for( std::list<M::CHalfEdge*>::iterator hiter = pL->halfedges().begin();
		  hiter != pL->halfedges().end(); hiter ++ )
	  {
		  M::CHalfEdge * pH = *hiter;

		  M::CEdge * pE = m_pHalfMesh->halfedgeEdge( pH );
		  M::CEdge * pWE = _match_ancestor_edge( pH );

		  if( pWE != NULL )
			  pE->length() =  pWE->length();
		  else
			  swapped_hes.push_back( pH );
	  }

	  //ph has been swapped
	  for( size_t j = 0; j < swapped_hes.size(); j ++ )
	  {
		  M::CHalfEdge * pH = swapped_hes[j];
		  _trace_swapped_edge( pH );
	  }
  }
};

/*!
 *	pH is a boundary halfedge on the half mesh, find the length of its
 *	ancestor edge on the double mesh
 */

template<typename M>
typename M::CEdge * CDynamicYamabeFlowPolyAnnulus<M>::_match_ancestor_edge( typename M::CHalfEdge * pH )
{
	M::CVertex * vs = m_pHalfMesh->halfedgeSource( pH );
	M::CVertex * vt = m_pHalfMesh->halfedgeTarget( pH );

	M::CVertex * ws = vs->ancestor();
	M::CVertex * wt = vt->ancestor();

	//the ancestor edge is on the double mesh
	for( M::VertexEdgeIterator veiter( ws); !veiter.end(); veiter ++ )
	{
		M::CEdge   *  e = *veiter;
		M::CVertex * w1 = m_pMesh->edgeVertex1( e ); 
		M::CVertex * w2 = m_pMesh->edgeVertex2( e );

		if( w1 == ws && w2 == wt ) return e;
		if( w1 == wt && w2 == ws ) return e;
	}

	return NULL;
};

/*!
 *	pH is a boundary halfedge on the half mesh, which has been swapped
 *  trace the symmetry polygon on the double mesh
 *
 */

template<typename M>
void CDynamicYamabeFlowPolyAnnulus<M>::_trace_swapped_edge( typename M::CHalfEdge * pH )
{
	M::CVertex * vs = m_pHalfMesh->halfedgeSource( pH );
	M::CVertex * vt = m_pHalfMesh->halfedgeTarget( pH );
	M::CEdge   * pe = m_pHalfMesh->halfedgeEdge( pH );

	for( std::set<M::CEdge*>::iterator eiter = m_cross_edges.begin(); eiter != m_cross_edges.end(); eiter ++ )
	{
		M::CEdge   * e = *eiter;
		M::CVertex * w = _corner_edge( e );
		if( w != vs && w != vt ) continue;
		double l = _trace( e );
		pe->length() = l;
	}
};




template<typename M>
double CDynamicYamabeFlowPolyAnnulus<M>::_trace( typename M::CEdge * head )
{
	std::vector<M::CHalfEdge*>  half_edges;
	
	M::CVertex * v_start = NULL;
	M::CVertex * v_end   = NULL;

	M::CEdge     * pe  = head;
	M::CHalfEdge * ph  = m_pMesh->edgeHalfedge( head, 0 );
	M::CVertex   * pv  = _third_vertex( ph );

	if( pv->dual() != pv )
	{
		ph = m_pMesh->edgeHalfedge( head, 1 );
		pv  = _third_vertex( ph );
		assert( pv->dual() == pv );
	}
	
	half_edges.push_back(ph);
	m_cross_edges.erase( head );

	v_start = pv;
		
	std::vector<M::CEdge*> swappable_edges;

	while( true )
	{

		M::CHalfEdge * sh = m_pMesh->halfedgeSym( ph );
		pv = _third_vertex( sh );

		if( pv->dual() == pv )
		{
			v_end = pv;
			_embed_symmetric_polygon( half_edges );
			CPoint2 d = v_end->huv() - v_start->huv();
			return d.norm();
		}
		
		M::CHalfEdge * nh = _parallel_transport_upper( ph );
		pe = m_pMesh->halfedgeEdge( nh );

		if( m_cross_edges.find( pe ) != m_cross_edges.end() )
		{	
			half_edges.push_back( nh );
			m_cross_edges.erase( pe );
			ph = nh;
			continue;
		}

		nh = _parallel_transport_lower( ph );
		pe = m_pMesh->halfedgeEdge( nh );

		if( m_cross_edges.find( pe ) != m_cross_edges.end() )
		{		
			half_edges.push_back( nh );
			m_cross_edges.erase( pe );
			ph = nh;
			continue;
		}

		assert(0);
	}
};


/*!
 *	embed the symmetric polygons along the symmetry axis
 */
template<typename M>
void CDynamicYamabeFlowPolyAnnulus<M>::_embed_symmetric_polygon( std::vector<typename M::CHalfEdge*> & half_edges )
{
	M::CVertex * v_start, *v_end;

	v_start = _third_vertex( half_edges[0]   );
	size_t n = half_edges.size();
	v_end   = _third_vertex( m_pMesh->halfedgeSym( half_edges[n-1] ) );

	std::queue<M::CFace*> faces;
	std::vector<M::CFace*> patch_faces;

	for( M::MeshVertexIterator viter( m_pMesh ) ; !viter.end(); viter ++ )
	{
		M::CVertex * pv = *viter;
		pv->touched() = false;
	}


	for( size_t i = 0; i < half_edges.size(); i ++ )
	{
		M::CHalfEdge * h = half_edges[i];
		faces.push( m_pMesh->halfedgeFace( h ) );
		patch_faces.push_back( m_pMesh->halfedgeFace( h ) );
		M::CHalfEdge * s = m_pMesh->halfedgeSym( h );
		faces.push( m_pMesh->halfedgeFace( s ) );
		patch_faces.push_back( m_pMesh->halfedgeFace( s ) );
	}	

	M::CFace * head = faces.front();
	faces.pop();
/*	
	for( FaceVertexIterator fviter( head ); !fviter.end(); fviter ++ )
	{
		V * pv = *fviter;
		std::cout << pv->id() << std::endl;
	}
*/
	_embed_first_face( head );

	while( !faces.empty() )
	{
		M::CFace * head = faces.front();
		faces.pop();
/*
		for( FaceVertexIterator fviter( head ); !fviter.end(); fviter ++ )
		{
			V * pv = *fviter;
			std::cout << pv->id() << std::endl;
		}
*/
		_embed_face( head );
	}

		//normalize the embedding
	std::set<M::CVertex*> patch_vertices;
		for( size_t i = 0; i < patch_faces.size(); i ++ )
		{
			M::CFace * f = patch_faces[i];

			for( M::FaceVertexIterator vfiter( f ); !vfiter.end() ; vfiter ++ )
			{
				M::CVertex * v = *vfiter;
				patch_vertices.insert( v );
			}
		}
		
		CPoint2 origin = v_start->huv();

		for( std::set<M::CVertex*>::iterator viter = patch_vertices.begin(); viter!= patch_vertices.end(); viter ++ )
		{
			M::CVertex * pv = *viter;
			pv->huv() = pv->huv() - origin;
		}

		CPoint2 d = v_end->huv() - v_start->huv();
		std::complex<double> r( d[0], d[1] );
		double theta = -std::arg(r);

		r  = std::complex<double>( cos(theta), sin(theta) );

		for( std::set<M::CVertex*>::iterator viter = patch_vertices.begin(); viter != patch_vertices.end(); viter ++ )
		{
			M::CVertex * pv = *viter;
			std::complex<double> z(pv->huv()[0], pv->huv()[1]);
			z = r * z;
			pv->huv() = CPoint2( z.real(), z.imag() );
		}
		//output embedding result
		//_debug_output( patch_faces );
}

//embed the first face
template<typename M>
void CDynamicYamabeFlowPolyAnnulus<M>::_embed_first_face( typename M::CFace * head )
{
	M::CHalfEdge * he[3];

	he[0] = m_pMesh->faceMostCcwHalfEdge( head  );
	he[1] = m_pMesh->faceNextCcwHalfEdge( he[0] );
	he[2] = m_pMesh->faceNextCcwHalfEdge( he[1] );

	std::vector<M::CVertex*> av;

	for( int i = 0; i < 3; i ++ )
	{
		av.push_back( m_pMesh->halfedgeTarget( he[(i+2)%3] ) );
	}

	M::CEdge * e = m_pMesh->halfedgeEdge( he[0] );
	av[0]->huv() = CPoint2(0,0);
	av[1]->huv() = CPoint2( e->length(), 0 );

	CPoint2 c1,c2;

	_circle_circle_intersection( CCircle( av[0]->huv(), m_pMesh->halfedgeEdge(he[2])->length()),
							     CCircle( av[1]->huv(), m_pMesh->halfedgeEdge(he[1])->length()),
								 c1,c2);

	if( cross( av[1]->huv() - av[0]->huv(), c1- av[0]->huv() ) > 0 )
	{
		av[2]->huv() = c1;
	}
	else
	{
		av[2]->huv() = c2;
	}


	for(int i = 0; i < 3; i ++ )
	{
		m_pMesh->halfedgeTarget( he[i] )->touched() = true;
	}

};

//vertex A, vertex B are known, vertex C is unknown

//embed the first face
template<typename M>
void CDynamicYamabeFlowPolyAnnulus<M>::_embed_face( typename M::CFace * head )
{
	std::vector<M::CVertex *> av;
	std::vector<M::CHalfEdge *> ah;

	for(M::FaceHalfedgeIterator fhiter(head); !fhiter.end(); fhiter ++ )
	{
		M::CHalfEdge * pH = *fhiter;
		ah.push_back( pH );
		M::CHalfEdge * pN = m_pMesh->halfedgeNext( pH );
		av.push_back( m_pMesh->halfedgeTarget( pN ) );
	}
/*
  for( FaceVertexIterator fviter( head ); !fviter.end(); fviter ++ )
  {
	  V * pV = *fviter;
	  av.push_back( pV );
  }
*/
	M::CVertex * A = NULL;
	M::CVertex * B = NULL;
	M::CVertex * C = NULL;

	M::CEdge * a = NULL;
	M::CEdge * b = NULL;
	M::CEdge * c = NULL;


  for( int i = 0; i < 3; i ++ )
  {
	  if( av[i]->touched() ) continue;
	  C = av[(i+0)%3]; c = m_pMesh->halfedgeEdge( ah[(i+0)%3] );
	  A = av[(i+1)%3]; a = m_pMesh->halfedgeEdge( ah[(i+1)%3] );
	  B = av[(i+2)%3]; b = m_pMesh->halfedgeEdge( ah[(i+2)%3] );
	  break;
  }
	
  if( C == NULL ) return;

	//radius of the first circle
	double r1 = b->length(); //m_pMesh->vertexEdge(A,C)->length();
	//radius of the second circle
	double r2 = a->length(); //m_pMesh->vertexEdge(B,C)->length();
	//center of the first circle
	CPoint2 c1 = A->huv();
	//center of the second circle
	CPoint2 c2 = B->huv();

	
	CPoint2 i1;
	CPoint2 i2;

	int result =_circle_circle_intersection( CCircle(c1,r1), CCircle(c2,r2), i1, i2);


	if( result )
	{
		if( cross( c2-c1, i1 - c1 ) > 0 )
			C->huv()  = i1;
		else
			C->huv()  = i2;

		C->touched() = true;
	}

	
	// This should never happen, once this happends,
	// we use crude approximation

	if( ! result ) 
	{
		assert(0);
	}
};



/*!
 *	inverse process of double covering
 */
/*
template<typename M>
void CDynamicYamabeFlowPolyAnnulus<M>::_debug_output( std::vector<F*>& patch_faces )
{

		static int idx = 0;
		std::string output("test_patch_");
		
		output += idx;
		output +=".m";

		std::fstream _os( output.c_str(), std::fstream::out );
		if( _os.fail() )
		{
			std::cerr<<"Error is opening file " << output << std::endl;
			return;
		}

		//normalize the embedding
		std::set<V*> patch_vertices;
		for( size_t i = 0; i < patch_faces.size(); i ++ )
		{
			F * f = patch_faces[i];

			for( FaceVertexIterator vfiter( f ); !vfiter.end() ; vfiter ++ )
			{
				V * v = *vfiter;
				patch_vertices.insert( v );
			}
		}

		for( std::set<V*>::iterator viter = patch_vertices.begin(); viter != patch_vertices.end(); viter ++ )
		{
			V * v = *viter;
			_os << "Vertex " << v->id() << " " << v->huv()[0] << " " << v->huv()[1]<< " 0" << std::endl;
		}
		for( size_t i = 0; i < patch_faces.size(); i ++ )
		{
			F * f = patch_faces[i];

			int j = 0; 
			int id[3];
			for( FaceVertexIterator vfiter( f ); !vfiter.end() ; vfiter ++ )
			{
				V * v = *vfiter;
				id[j++] = v->id();
			}
			_os << "Face " << f->id() << " " << id[0] << " " << id[1]<< " " << id[2] << std::endl;
		}
		_os.close();
};
*/

template<typename M>
void CDynamicYamabeFlowPolyAnnulus<M>::_calculate_metric_circle_domain_on_cylinder()
{
	  for( int i = 0; i < 32; i ++ )
	  {
		  _set_target_curvature();
		  double error_threshold = 1e-10;
		  double step_length = 1.0;
		  _Dynamic_Newton( error_threshold, step_length );
	  }
};

} //namespace RicciFlow
} //namespace MeshLib

#endif  _DYNAMIC_YAMABE_FLOW_POLY_ANNULUS_H_