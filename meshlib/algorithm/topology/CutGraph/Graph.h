#ifndef _GRAPH_H_
#define _GRAPH_H_
#include <map>
#include <vector>

#include "Mesh/BaseMesh.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "mesh/boundary.h"
#include "Parser/parser.h"

namespace MeshLib
{
namespace Topology
{

template<typename M>
class	CGraph
{
public:
	CGraph( M * pMesh ){ m_pMesh=pMesh; m_root=NULL;};
	~CGraph()
	{
		for( std::list<std::list<M::CHalfEdge*>*>::iterator liter = m_loops.begin(); liter != m_loops.end(); liter ++ )
	  {
		  std::list<M::CHalfEdge*>* pL = *liter;
		delete pL;
	  }	
	};

	void  _insert( typename M::CEdge * pEdge );
  void  _breadth_first_search();
  void  _locate_loops();

  std::list<std::list<typename M::CHalfEdge*>*>  & loops() 
  {
    return m_loops;
  };
protected:
  M * m_pMesh;
  std::set<typename M::CEdge*  >   m_edges;
  std::set<typename M::CVertex*>   m_vertices;
  std::list<typename M::CEdge* >   m_left_edges;
  std::list<std::list<typename M::CHalfEdge*>*>  m_loops;
  typename M::CVertex * m_root;

  void _trace_to_root( typename M::CVertex * leaf, typename M::CVertex * root, std::vector<typename M::CVertex*> & path );
};

template<typename M>
void CGraph<M>::_insert( typename M::CEdge * pEdge )
{
    m_edges.insert( pEdge );
	M::CVertex * v1 = m_pMesh->edgeVertex1( pEdge );
	M::CVertex * v2 = m_pMesh->edgeVertex2( pEdge );
    m_vertices.insert( v1 );
    m_vertices.insert( v2 );
};

template<typename M>
void CGraph<M>::_breadth_first_search()
{
	for( std::set<M::CVertex*>::iterator viter = m_vertices.begin(); viter != m_vertices.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
        v->touched() = false;
        v->father() = NULL;
        m_root = v;
  }
  
	std::queue<M::CVertex*> vqueue;

  vqueue.push( m_root );
  m_root->touched() = true;
  m_root->father()  = NULL;

  while( !vqueue.empty() )
  {
	  M::CVertex * head = vqueue.front();
      vqueue.pop();

	  for( M::VertexEdgeIterator veiter( head );  !veiter.end(); ++veiter )
      {
		  M::CEdge * e    = *veiter;
		  std::set<M::CEdge*>::iterator eiter = m_edges.find( e );
          if( eiter == m_edges.end() ) continue;
		  M::CVertex * w = (m_pMesh->edgeVertex1(e ) != head )? m_pMesh->edgeVertex1(e): m_pMesh->edgeVertex2(e);
          if( w->touched() ) continue;
          w->father() = head;
          vqueue.push( w );
          w->touched() = true;
      }
  }
  std::set<M::CEdge*> edges = m_edges;

  for( std::set<M::CVertex*>::iterator viter = m_vertices.begin(); viter != m_vertices.end(); viter ++ )
    {
		M::CVertex * v = *viter;
		M::CVertex * w = v->father();
        if( w == NULL ) continue;
		M::CEdge * e = m_pMesh->vertexEdge( v,w);
		std::set<M::CEdge*>::iterator epos = edges.find( e );
        assert( epos != edges.end() );
        edges.erase( epos );
    }

    printf("Left %d edges\n", edges.size());
    
	for( std::set<M::CEdge*>::iterator eiter = edges.begin(); eiter != edges.end(); eiter ++ )
    {
		M::CEdge * e = *eiter;
		m_left_edges.push_back( e );
    }
  
};

template<typename M>
void CGraph<M>::_trace_to_root( typename M::CVertex * leaf, typename M::CVertex * root, std::vector<typename M::CVertex *> & path )
{
    path.push_back( leaf );
	M::CVertex * father = leaf->father();
    while( father != NULL )
    {
        path.push_back( father );
        father = father->father();
    }
};

//compute the loops in the graph, 

template<typename M>
void CGraph<M>::_locate_loops()
{
	for( std::list<M::CEdge*>::iterator eiter = m_left_edges.begin(); eiter != m_left_edges.end(); eiter ++ )
  {
	  M::CEdge * e = *eiter;
	  M::CVertex * v1 = m_pMesh->edgeVertex1( e );
	  M::CVertex * v2 = m_pMesh->edgeVertex2( e );
  
	  std::vector<M::CVertex*> trace1;
	  std::vector<M::CVertex*> trace2;

      _trace_to_root( v1, m_root, trace1 );
      _trace_to_root( v2, m_root, trace2 );

	  int s1 = trace1.size()-1;
	  int s2 = trace2.size()-1;
	  while( trace1[s1] == trace2[s2] )
	  {
		  s1 --; s2 --;
	  }
	  s1++;s2++;

	  std::vector<M::CVertex*> path1;
	  std::vector<M::CVertex*> path2;

	  for( int j = 0 ; j < s1 + 1 ; j ++ )
	  {
		  path1.push_back( trace1[j] );
	  }

	  for( int j = 0 ; j < s2 + 1; j ++ )
	  {
		  path2.push_back( trace2[j] );
	  }

      //merge two paths to form a loop

	  std::list<M::CHalfEdge*> * pLoop = new std::list<M::CHalfEdge*>;
      assert( pLoop );
      m_loops.push_back( pLoop );

      for( int i = path1.size()-1; i > 0; i-- )
      {
		  M::CVertex * s = path1[i];
		  M::CVertex * t = path1[i-1];
		  M::CEdge *     e = m_pMesh->vertexEdge( s,t );
		  M::CHalfEdge * he = e->halfedge(0);
          if( he->target() != t )
          {
              he = e->halfedge(1);
          }
          assert( he->target() == t && he->source() == s );
          pLoop->push_back( he );
      }

	  M::CHalfEdge * he = e->halfedge(0);
      if( he->target () != v2 )
      {
        he = e->halfedge(1);
      }

      pLoop->push_back( he );

      for( size_t i = 0 ; i < path2.size()-1; i++ )
      {
		  M::CVertex * s = path2[i];
		  M::CVertex * t = path2[i+1];
		  M::CEdge   * e = m_pMesh->vertexEdge( s,t );
		  M::CHalfEdge * he = e->halfedge(0);
          if( he->target() != t )
          {
              he = e->halfedge(1);
          }
          assert( he->target() == t && he->source() == s );
          pLoop->push_back( he );
      }

      printf("Loop\n");
	  for( std::list<M::CHalfEdge*>::iterator hiter = pLoop->begin(); hiter != pLoop->end(); hiter++ )
      {
		  M::CHalfEdge * he = *hiter;
          printf("%d -> %d ", he->source()->id(), he->target()->id() );
      }
    printf("\n\n");
  };
};

//typedef CGraph<CCutGraphVertex, CCutGraphEdge, CCutGraphFace, CHalfEdge> Graph;

} //namespace Topology
} //namespace MeshLib

#endif  //_GRAPH_H_