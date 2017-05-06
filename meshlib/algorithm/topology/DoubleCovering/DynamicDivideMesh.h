#ifndef  _DYNAMIC_DIVIDE_MESH_H_
#define  _DYNAMIC_DIVIDE_MESH_H_

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
namespace Topology
{
/*!
 *	Compare vertex father
 */
template<typename V>
class CompareVertexFather {
public:
    bool operator()(V* & v1, V* & v2)
    {
       if (v1->father() > v2->father() ) return true;
       if (v1->father() < v2->father() ) return false;
       return true;
    }
};


/*! \brief CDynamicDivideMesh class
*   
*  Divide the doubled mesh to two halves, split all the edges cross the symmetry axis
*
*/
template<typename V, typename E, typename F, typename H>
class CDynamicDivideMesh : public CDynamicMesh<V,E,F,H>
{
public:
	typedef V CVertex;
	typedef E CEdge;
	typedef F CFace;
	typedef H CHalfEdge;

	typedef MeshVertexIterator<V,E,F,H>        MeshVertexIterator;
	typedef MeshFaceIterator<V,E,F,H>          MeshFaceIterator;
	typedef MeshEdgeIterator<V,E,F,H>          MeshEdgeIterator;
	typedef FaceVertexIterator<V,E,F,H>        FaceVertexIterator;
	typedef FaceHalfedgeIterator<V,E,F,H>	   FaceHalfedgeIterator;
	typedef VertexEdgeIterator<V,E,F,H>		   VertexEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H>	   VertexVertexIterator;
	typedef FaceEdgeIterator<V,E,F,H>		   FaceEdgeIterator;
	

	/*!
	 *	Divide the mesh
	 */
	void Divide( const char * output );

	void _swap_edge( E * e );

private:
	H * _parallel_transport_upper( H * h )
	{
		return halfedgeNext( halfedgeSym( halfedgePrev( halfedgeSym(h) ) ) );
	}
	H * _parallel_transport_lower( H * h )
	{
		return halfedgePrev( halfedgeSym( halfedgeNext( halfedgeSym(h) ) ) );
	}
	V * _third_vertex( H * h )
	{
		return halfedgeTarget( halfedgeNext(h) );
	}
private:
	//embedding a sequence of ajdacent faces

	//embed the first face
	void _embed_first_face( F * head );
	//embed other face
	void _embed_face( F * head );	
	/*!
	 *	embed a symmetric polygon
	 */
	void _embed_symmetric_polygon( std::vector<H*> & half_edges );

private:
	//e is a cross edge, and e is against a self-dual vertex
	bool _corner_edge( E* e );
	/*!
	 *	match vertex with its dual
	 */
	void _match_dual_vertex();

private:
	
	//trace from a corner edge to get a symmetric polygon
	void _trace( std::set<E*> & crosses, E * e );
	//trace from a corner edge to get a symmetric polygon
	void _remesh( std::set<E*> & crosses, E * e );

	/*!
	 *	output the faces for debug purposes
	 */
	void _debug_output( std::vector<F*> & faces );
	/*!
	 *	combinatorial structure surgery, to make the symmetry axis along with the edges
	 */
	void _surgery( std::vector<H*>& half_edges );
	/*!
	 *	split the doubled mesh to 2 halves, output one half
	 */
	void _split( const char * output );

};

/*!
 *	embed the symmetric polygons along the symmetry axis
 */
template<typename V, typename E, typename F, typename H>
void CDynamicDivideMesh<V,E,F,H>::_embed_symmetric_polygon( std::vector<H*> & half_edges )
{
	V * v_start, *v_end;

	v_start = _third_vertex( half_edges[0]   );
	size_t n = half_edges.size();
	v_end   = _third_vertex( halfedgeSym( half_edges[n-1] ) );

	std::queue<F*> faces;
	std::vector<F*> patch_faces;

	for( MeshVertexIterator viter( this ) ; !viter.end(); viter ++ )
	{
		V * pv = *viter;
		pv->touched() = false;
	}


	for( size_t i = 0; i < half_edges.size(); i ++ )
	{
		H * h = half_edges[i];
		faces.push( halfedgeFace( h ) );
		patch_faces.push_back( halfedgeFace( h ) );
		H * s = halfedgeSym( h );
		faces.push( halfedgeFace( s ) );
		patch_faces.push_back( halfedgeFace( s ) );
	}	

	F * head = faces.front();
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
		F * head = faces.front();
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
		
		CPoint2 origin = v_start->huv();

		for( std::set<V*>::iterator viter = patch_vertices.begin(); viter!= patch_vertices.end(); viter ++ )
		{
			V * pv = *viter;
			pv->huv() = pv->huv() - origin;
		}

		CPoint2 d = v_end->huv() - v_start->huv();
		std::complex<double> r( d[0], d[1] );
		double theta = -std::arg(r);

		r  = std::complex<double>( cos(theta), sin(theta) );

		for( std::set<V*>::iterator viter = patch_vertices.begin(); viter != patch_vertices.end(); viter ++ )
		{
			V * pv = *viter;
			std::complex<double> z(pv->huv()[0], pv->huv()[1]);
			z = r * z;
			pv->huv() = CPoint2( z.real(), z.imag() );
		}
		//output embedding result
		_debug_output( patch_faces );
}

/*!
 *	inverse process of double covering
 */
template<typename V, typename E, typename F, typename H>
void CDynamicDivideMesh<V,E,F,H>::_debug_output( std::vector<F*>& patch_faces )
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

/*!
 *	two vertices are dual to each other, if they share the same father,
 *	Assume the father field has been set
 *  The vertex dual field will be set
 */
template<typename V, typename E, typename F, typename H>
void CDynamicDivideMesh<V,E,F,H>::_match_dual_vertex()
{
	std::priority_queue<V*, std::vector<V*>, CompareVertexFather<V> > pq;
	for( MeshVertexIterator viter( this ); !viter.end(); viter ++ )
	{
		V * pV = *viter;
		pq.push( pV );
	}
	std::vector<V*> vs;
	while (! pq.empty()) {
		V * pV = pq.top();
		vs.push_back( pV );
		pq.pop();
	}
	for( size_t i = 0; i < vs.size(); i ++ )
	{
		V * pV = vs[i];
		V * pW = vs[(i+1)%vs.size()];

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

/*!
 *	inverse process of double covering
 */
template<typename V, typename E, typename F, typename H>
void CDynamicDivideMesh<V,E,F,H>::Divide( const char * output )
{

	_match_dual_vertex();

	//edges crosses the boader, including vertical edges for diamonds, and rectangular shapes
	std::set<E*> crosses;

	for( MeshEdgeIterator eiter( this ); !eiter.end(); eiter ++ )
	{
		E * e = *eiter;
		V * v1 = edgeVertex1( e );
		V * v2 = edgeVertex2( e );
		if( v1 == v2 ) continue;

		if( v1->dual() == v2 && v2->dual() == v1 )
			crosses.insert( e );
	}

	while( !crosses.empty() )
	{
		E * head = NULL;
		for( std::set<E*>::iterator eiter = crosses.begin(); eiter != crosses.end(); eiter ++ )
		{
			head = *eiter;
			if( _corner_edge( head ) ) break;
		}
		//trace the corner edge to get a symmetric polygon
		if( head != NULL )
		{
			_trace( crosses, head );
			std::cout << crosses.size() << std::endl;
		}
		else
		{
			assert(0);
			break;
		}
	}


	//edges crosses the boader, including vertical edges for diamonds, and rectangular shapes

	for( MeshEdgeIterator eiter( this ); !eiter.end(); eiter ++ )
	{
		E * e = *eiter;

		V * v1 = edgeVertex1( e );
		V * v2 = edgeVertex2( e );

		if( v1 == v2 ) continue;

		if( v1->dual() == v2 && v2->dual() == v1 )
			crosses.insert( e );
	}

	while( !crosses.empty() )
	{
		E * head = NULL;
		for( std::set<E*>::iterator eiter = crosses.begin(); eiter != crosses.end(); eiter ++ )
		{
			head = *eiter;
			if( _corner_edge( head ) ) break;
		}
		//trace the corner edge to get a symmetric polygon
		if( head != NULL )
		{
			_remesh( crosses, head );
			std::cout << crosses.size() << std::endl;
		}
		else
		{
			assert(0);
			break;
		}
	}
	
	_split( output );
};

template<typename V, typename E, typename F, typename H>
void CDynamicDivideMesh<V,E,F,H>::_swap_edge( E * e )
{
		H *    lh[3];
		double ll[3];
		double la[3];

		lh[0] = edgeHalfedge( e, 0 );

		for( int i = 0; i < 2; i ++ )
		{
			lh[i+1] = faceNextCcwHalfEdge( lh[i] );
		}

		for( int i = 0; i < 3; i ++ )
		{
			E * e = halfedgeEdge( lh[i] );
			ll[i] = e->length();
		}

		for( int i = 0; i < 3; i ++ )
		{
			la[i] = ( ll[(i+0)%3]*ll[(i+0)%3] + ll[(i+1)%3]* ll[(i+1)%3] - ll[(i+2)%3]*ll[(i+2)%3] )/(2 * ll[(i+0)%3] * ll[(i+1)%3]);
		}

		H *    rh[3];
		double rl[3];
		double ra[3];

		rh[0] = edgeHalfedge( e, 1 );

		for( int i = 0; i < 2; i ++ )
		{
			rh[i+1] = faceNextCcwHalfEdge( rh[i] );
		}

		for( int i = 0; i < 3; i ++ )
		{
			E * e = halfedgeEdge( rh[i] );
			rl[i] = e->length();
		}

		for( int i = 0; i < 3; i ++ )
		{
			ra[i] = ( rl[(i+0)%3]*rl[(i+0)%3] + rl[(i+1)%3]* rl[(i+1)%3] - rl[(i+2)%3]*rl[(i+2)%3] )/(2 * rl[(i+0)%3] * rl[(i+1)%3]);
		}

		double diagonal = sqrt( ll[1]*ll[1] + rl[2] * rl[2] - 2 * ll[1] * rl[2] * cos( la[0] + ra[2] ) );
		std::cerr << "Edge Swap" << std::endl;
		swapEdge( e );
		//e->length() = diagonal;

};

/*!
 *	an edge is a corner edge if its 3rd vertex is self-dual
 */
template<typename V, typename E, typename F, typename H>
bool CDynamicDivideMesh<V,E,F,H>::_corner_edge( E * e )
{
	H * h = edgeHalfedge( e, 0 );
	V * v = _third_vertex( h );
	if( v->dual() == v ) return true;

	h = edgeHalfedge( e, 1 );
	v = _third_vertex( h );
	if( v->dual() == v ) return true;
	
	return false;
};

/*!
 *	from a corner edge trace a sequence of cross edges, which are parallel to each other
 *
 *	an edge is a cross edge if its end verticies are dual to each other
 *
 */

template<typename V, typename E, typename F, typename H>
void CDynamicDivideMesh<V,E,F,H>::_trace( std::set<E*> & crosses, E * head )
{
	std::vector<E*>  cross_edges;
	std::vector<H*>  half_edges;
	
	V * v_start = NULL;
	V * v_end   = NULL;

	E * pe  = head;
	H * ph  = edgeHalfedge( head, 0 );
	V * pv  = _third_vertex( ph );

	if( pv->dual() != pv )
	{
		ph = edgeHalfedge( head, 1 );
		pv  = _third_vertex( ph );
		assert( pv->dual() == pv );
	}
	
	cross_edges.push_back( pe );
	half_edges.push_back(ph);
	crosses.erase( head );

	v_start = pv;
		
	std::vector<E*> swappable_edges;

	while( true )
	{

		H * sh = halfedgeSym( ph );
		pv = _third_vertex( sh );

		if( pv->dual() == pv )
		{
			v_end = pv;

			//_embed_symmetric_polygon( half_edges );


			for( size_t i = 0; i < swappable_edges.size(); i ++ )
			{
				E * e = swappable_edges[i];
				_swap_edge( e );
			}

			return;
		}
		
		H * nh = _parallel_transport_upper( ph );
		pe = halfedgeEdge( nh );

		if( crosses.find( pe ) != crosses.end() )
		{	
			cross_edges.push_back( pe );
			half_edges.push_back( nh );
			crosses.erase( pe );
			ph = nh;
			//swap these edges
			swappable_edges.push_back( halfedgeEdge( halfedgePrev(nh) ) );
			continue;
		}

		nh = _parallel_transport_lower( ph );
		pe = halfedgeEdge( nh );

		if( crosses.find( pe ) != crosses.end() )
		{		
			cross_edges.push_back( pe );
			half_edges.push_back( nh );
			crosses.erase( pe );
			ph = nh;
			continue;
		}

		assert(0);
	}

};


/*!
 *	from a corner edge trace a sequence of cross edges, which are parallel to each other
 *
 *	an edge is a cross edge if its end verticies are dual to each other
 *
 */

template<typename V, typename E, typename F, typename H>
void CDynamicDivideMesh<V,E,F,H>::_remesh( std::set<E*> & crosses, E * head )
{
	std::vector<E*>  cross_edges;
	std::vector<H*>  half_edges;
	//std::vector<E*>  diagonals;

	V * v_start = NULL;
	V * v_end   = NULL;

	E * pe  = head;
	H * ph  = edgeHalfedge( head, 0 );
	V * pv  = _third_vertex( ph );

	if( pv->dual() != pv )
	{
		ph = edgeHalfedge( head, 1 );
		pv  = _third_vertex( ph );
		assert( pv->dual() == pv );
	}
	
	cross_edges.push_back( pe );
	half_edges.push_back(ph);
	crosses.erase( head );

	v_start = pv;
		

	while( true )
	{

		H * sh = halfedgeSym( ph );
		pv = _third_vertex( sh );

		if( pv->dual() == pv )
		{
			v_end = pv;

			_embed_symmetric_polygon( half_edges );
			_surgery( half_edges );	
			return;
		}
		
		H * nh = _parallel_transport_upper( ph );
		nh = _parallel_transport_lower( ph );
		pe = halfedgeEdge( nh );

		if( crosses.find( pe ) != crosses.end() )
		{		
			cross_edges.push_back( pe );
			half_edges.push_back( nh );
			crosses.erase( pe );
			ph = nh;
			continue;
		}

		assert(0);
	}

};


//embed the first face
template<typename V, typename E, typename F, typename H>
void CDynamicDivideMesh<V,E,F,H>::_embed_first_face( F * head )
{
	H * he[3];

	he[0] = faceMostCcwHalfEdge( head  );
	he[1] = faceNextCcwHalfEdge( he[0] );
	he[2] = faceNextCcwHalfEdge( he[1] );

	std::vector<V*> av;

	for( int i = 0; i < 3; i ++ )
	{
		av.push_back( halfedgeTarget( he[(i+2)%3] ) );
	}

	E * e = halfedgeEdge( he[0] );
	av[0]->huv() = CPoint2(0,0);
	av[1]->huv() = CPoint2( e->length(), 0 );

	CPoint2 c1,c2;

	_circle_circle_intersection( CCircle( av[0]->huv(), halfedgeEdge(he[2])->length()),
							     CCircle( av[1]->huv(), halfedgeEdge(he[1])->length()),
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
		halfedgeTarget( he[i] )->touched() = true;
	}

};

//vertex A, vertex B are known, vertex C is unknown

//embed the first face
template<typename V, typename E, typename F, typename H>
void CDynamicDivideMesh<V,E,F,H>::_embed_face( F * head )
{
	std::vector<V*> av;
	std::vector<H*> ah;

	for(FaceHalfedgeIterator fhiter(head); !fhiter.end(); fhiter ++ )
	{
		H * pH = *fhiter;
		ah.push_back( pH );
		H * pN = halfedgeNext( pH );
		av.push_back( halfedgeTarget( pN ) );
	}
/*
  for( FaceVertexIterator fviter( head ); !fviter.end(); fviter ++ )
  {
	  V * pV = *fviter;
	  av.push_back( pV );
  }
*/
  V * A = NULL;
  V * B = NULL;
  V * C = NULL;

  E * a = NULL;
  E * b = NULL;
  E * c = NULL;


  for( int i = 0; i < 3; i ++ )
  {
	  if( av[i]->touched() ) continue;
	  C = av[(i+0)%3]; c = halfedgeEdge( ah[(i+0)%3] );
	  A = av[(i+1)%3]; a = halfedgeEdge( ah[(i+1)%3] );
	  B = av[(i+2)%3]; b = halfedgeEdge( ah[(i+2)%3] );
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



//embed the first face
template<typename V, typename E, typename F, typename H>
void CDynamicDivideMesh<V,E,F,H>::_surgery( std::vector<H*>& half_edges )
{

	V * v_start, *v_end;
	v_start = _third_vertex( half_edges[0]   );
	size_t n = half_edges.size();
	v_end   = _third_vertex( halfedgeSym( half_edges[n-1] ) );

	std::vector<E*> cross_edges;
	std::vector<E*> diagonals;

	for( size_t i = 0; i < half_edges.size(); i ++ )
	{
		cross_edges.push_back( halfedgeEdge( half_edges[i]) );
	}

	for( size_t i = 0; i < half_edges.size()-1; i ++ )
	{
		H * h = half_edges[i];
		E * e = halfedgeEdge( halfedgeNext( halfedgeSym( h ) ) );
		diagonals.push_back( e );
	}	

	m_vertex_id = 0;
	for( MeshVertexIterator viter( this ); !viter.end(); viter ++ )
	{
		V * pv = *viter;
		m_vertex_id = (m_vertex_id > pv->id() )? m_vertex_id: pv->id();
	}

		m_edge_id = 0;
		for( MeshEdgeIterator eiter( this ); !eiter.end(); eiter ++ )
		{
			E * pe = *eiter;
			m_edge_id = (m_edge_id > pe->id() )? m_edge_id: pe->id();
		}

		m_face_id = 0;
		for( MeshFaceIterator fiter( this ); !fiter.end(); fiter ++ )
		{
			F * pf = *fiter;
			m_face_id = (m_face_id > pf->id() )? m_face_id: pf->id();
		}

			
		std::vector<V*> new_verts;
		for( size_t i = 0; i < cross_edges.size(); i ++ )
		{
			E * pe = cross_edges[i];
			V * v1 = edgeVertex1( pe );
			V * v2 = edgeVertex2( pe );

			V * pv = splitEdge( pe );
			pv->huv() = (v1->huv() + v2->huv())/2.0;

			double rho = (v1->huv()[0] - v_start->huv()[0])/(v_end->huv()[0] - v_start->huv()[0]);

			pv->point() = v_start->point()*(1-rho) + v_end->point() * rho;
			pv->rgb()   = v_start->rgb()*(1-rho) + v_end->rgb() * rho;

			pv->father() = pv->id();
			pv->dual() = pv;
			new_verts.push_back( pv );
		}


		for( size_t i = 0; i < diagonals.size(); i ++ )
		{
			E * pe = diagonals[i];
			_swap_edge( pe );
		}

		for( size_t i = 0; i < new_verts.size(); i ++ )
		{
			V * pv = new_verts[i];
			for( VertexEdgeIterator veiter( pv ); !veiter.end(); veiter ++ )
			{
				E * e = *veiter;
				V * v1 = edgeVertex1( e );
				V * v2 = edgeVertex2( e );
				CPoint2 d = v1->huv() - v2->huv();
				e->length() = d.norm();
			}
		}
};


//embed the first face
template<typename V, typename E, typename F, typename H>
void CDynamicDivideMesh<V,E,F,H>::_split( const char * output )
{
	//locate the symmetry axis 
	for( MeshEdgeIterator eiter(this ); !eiter.end(); eiter ++ )
	{
		E * e = *eiter;
		e->sharp() = false;

		V * v1 = edgeVertex1( e );
		V * v2 = edgeVertex2( e );
		if( v1->dual() == v1 && v2->dual()== v2 )
		{
			e->sharp() = true;
		}
	}

	for( MeshFaceIterator fiter( this ); !fiter.end(); fiter ++ )
	{
		F * pf = *fiter;
		pf->touched() = false;
	}

	std::queue<F*> fqueue;

	for( MeshFaceIterator fiter( this ); !fiter.end(); fiter ++ )
	{
		F * pf = *fiter;
		fqueue.push( pf );
		pf->touched() = true;
		break;
	}

	std::set<F*> sfaces;

	while( !fqueue.empty() )
	{
		F * pf = fqueue.front();
		fqueue.pop();
		sfaces.insert( pf );
		//std::cout << pf->id() << std::endl;

		for( FaceHalfedgeIterator fhiter( pf ); !fhiter.end(); fhiter ++ )
		{
			H * h = *fhiter;
			E * e = halfedgeEdge( h );
			
			if( e->sharp() ) continue;

			H * s  = halfedgeSym(h);
			F * wf = halfedgeFace( s );
			//std::cout << wf->id() << std::endl;

			//if( sfaces.find( wf ) == sfaces.end() )
			if( !wf->touched() )
			{
				fqueue.push( wf );
				wf->touched() = true;
			}
		}
	}
/*
	for( MeshVertexIterator viter( this ); !viter.end(); viter ++ )
	{
		V * pv = *viter;
		pv->touched() = false;
	}

	std::queue<V*> vqueue;

	for( MeshVertexIterator viter( this ); !viter.end(); viter ++ )
	{
		V * pv = *viter;
		if( pv->dual() == pv )
		{
			pv->touched() = true;
			//vqueue.push( pv );
		}
	}

	for( MeshVertexIterator viter( this ); !viter.end(); viter ++ )
	{
		V * pv = *viter;
		if( ! pv->touched() )
		{
			vqueue.push(pv);
			pv->touched() = true;
			break;
		}
	}

	while( !vqueue.empty() )
	{
		V * pv = vqueue.front();
		vqueue.pop();

		for( VertexVertexIterator vviter(pv); !vviter.end(); vviter ++ )
		{
			V * pw = *vviter;
			if( pw->touched() ) continue;
			pw->touched() = true;
			vqueue.push( pw );
		}
	}
*/
	std::set<E*> edges;
/*
	std::set<F*> faces;

	for( MeshFaceIterator fiter(this); !fiter.end(); fiter ++ )
	{
		F * pf = *fiter;

		bool touched = true;
		for( FaceVertexIterator fviter( pf ); !fviter.end(); fviter ++ )
		{
			V * pv = *fviter;
			if( !pv->touched() ) touched = false;
		}
		if( !touched ) continue;
		faces.insert( pf );
	}
*/
	for( std::set<F*>::iterator fiter = sfaces.begin(); fiter != sfaces.end(); fiter ++ )
	{
		F * pf = *fiter;
		for( FaceEdgeIterator feiter( pf ); !feiter.end(); feiter ++ )
		{
			E * pe = *feiter;
			edges.insert( pe );
		}
	}

	std::set<V*> verts;
	for( std::set<F*>::iterator fiter = sfaces.begin(); fiter != sfaces.end(); fiter ++ )
	{
		F * pf = *fiter;
		for( FaceVertexIterator fviter( pf ); !fviter.end(); fviter ++ )
		{
			V * pv = *fviter;
			verts.insert( pv );
		}
	}





	//write traits to string
	for( MeshVertexIterator viter( this ); !viter.end(); viter ++ )
	{
		V * pV = *viter;
		pV->_to_string();
	}
	for( MeshEdgeIterator eiter(this); !eiter.end(); eiter ++ )
	{
		E * pE = *eiter;
		pE->_to_string();
	}
	for( MeshFaceIterator fiter( this ); !fiter.end(); fiter ++ )
	{
		F * pF = *fiter;
		pF->_to_string();
	}
	for( MeshFaceIterator fiter(this ); !fiter.end(); fiter ++ )
	{
		F * pF = *fiter;
		H * pH  = faceMostCcwHalfEdge( pF );
		do{
			pH->_to_string();
			pH = faceNextCcwHalfEdge( pH );
		}while( pH != faceMostCcwHalfEdge(pF ) );
	}

	std::fstream _os( output, std::fstream::out );
	if( _os.fail() )
	{
		fprintf(stderr,"Error is opening file %s\n", output );
		return;
	}

	for( std::set<V*>::iterator viter= verts.begin(); viter != verts.end(); viter ++ )
	{
		V * v = *viter;

		_os << "Vertex " << v->id();
		for( int i = 0; i < 3; i ++ )
		{
			_os << " " << v->point()[i];
		}
		if( v->string().size() > 0 )
		{
			_os << " " <<"{"<< v->string() << "}";
		}
		_os << std::endl;
	}

	int eidx = 1;

	for( std::set<E*>::iterator eiter = edges.begin(); eiter != edges.end(); eiter ++ )
	{
		E * e = *eiter;
		e->index() = eidx ++;
	}

	for( std::set<E*>::iterator eiter = edges.begin(); eiter != edges.end(); eiter ++ )
	{
		E * e = *eiter;
		//_os << "Edge "<<  e->id() << " " << edgeVertex1(e)->id() <<" " << edgeVertex2(e)->id() << " ";
		_os << "Edge "<<  e->index() << " " << edgeVertex1(e)->id() <<" " << edgeVertex2(e)->id() << " ";
		if( e->string().size() > 0 )
		{
			_os << "{" << e->string() << "}";
		}
		_os << std::endl;
	}
	for( std::set<CFace*>::iterator fiter = sfaces.begin(); fiter != sfaces.end(); fiter ++ )
	{
		F * f = *fiter;
		_os << "Face " << f->id();
		H * he = faceHalfedge( f );
		do{
			E * e = halfedgeEdge(he);
			if ( he == edgeHalfedge(e, 0) )		// pos
				//_os << " +" <<  e->id();
				_os << " +" <<  e->index();
			else							// neg
				//_os << " -" <<  e->id();
				_os << " -" <<  e->index();
			he = halfedgeNext( he );
		}while( he != faceHalfedge( f ) );
		if( f->string().size() > 0 )
		{
			_os << " " << "{"<< f->string() << "}";
		}
		_os << std::endl;
	}
	for( std::set<F*>::iterator fiter = sfaces.begin(); fiter != sfaces.end(); fiter ++  )
	{
		F * f = *fiter;
		H * he = faceHalfedge( f );
		do{
  			if( he->string().size() > 0 )
			  {
				  _os << "Corner "<< he->vertex()->id() << " " << f->id() << " ";
				  _os << "{" << he->string() << "}" << std::endl;
			  }
			  he = halfedgeNext( he );
		}while( he != faceHalfedge(f) );
	}

	_os.close();

};

typedef CDynamicDivideMesh<CDynamicDoubleCoveringVertex, CDynamicDoubleCoveringEdge, CDynamicDoubleCoveringFace, CDynamicDoubleCoveringHalfEdge> CDDDMesh;

} //namespace Topology
} //namespace MeshLib

#endif