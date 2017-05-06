/*!
*      \file DynamicShortestPath.h
*      \brief Algorithm for shortest paths
*	   \author David Gu
*      \date Document 10/12/2013
*
*		Compute the shortest path from one inner boundary loop to the exterior boundary loop. The input is a mesh with 
*       multiple boundary components; the output are meshes with shortest paths, labeled as sharp edges,  between boundary 
*       k to boundary 0 (exterior boundary).
*/


#ifndef _DYNAMIC_SHORTEST_PATH_H_
#define _DYNAMIC_SHORTEST_PATH_H_

#include  <math.h>
#include <queue>
#include "Mesh/boundary.h"
#include "Mesh/iterators.h"
#include "DynamicShortestPathMesh.h"

namespace MeshLib
{
/*!
 * \brief CDynamicShortestPath class
 * 
 * Compute the shortest path between one interior boundary component and the exterior boundary component
 * 
 */
template <typename M>
  class CDynamicShortestPath
  {
  public:
    /*!
	 * CDynamicShortestPath default constructor
	 * \param pMesh the input mesh
	 */
    CDynamicShortestPath( M * pMesh );
	/*! CShortestPath destructor
	 *
	 */
    ~CDynamicShortestPath();
	/*!	Compute the shortest path from boundary loop k to boundary loop 0 (exterior boundary loop ),
	 *  the path is labeled as sharp edges and output to "output_name_k.m"
	 *  \param prefix the prefix of output mesh name
	 */
    void _cut( const char * prefix );
  
  protected:
    /*! Pointer to the input mesh
	 */
    M * m_pMesh;
	/*!	Boundary loops of the input mesh
	 */
	typename M::CBoundary m_boundary;
	/*! Edges on the shortest path 
	 */
	std::list<typename M::CEdge*> m_cuts;
	/*!	Compute the shortest path between source loop to the target loop
	 * \param source the starting boundary loop
	 * \param target the ending boundary loop
	 */
	void _trace(typename M::CLoop * source, typename M::CLoop * target);
  };


/*! CDynamicShortestPath constructor
 *
 *  \param pMesh the input mesh
 */
template <typename M>
CDynamicShortestPath<M>::CDynamicShortestPath( M * pMesh ): m_pMesh( pMesh ), m_boundary( m_pMesh )
{
};

/*! CDynamicShortestPath destructor
 *
 */

template <typename M>
CDynamicShortestPath<M>::~CDynamicShortestPath()
{
};

/*!	Compute the shortest paths between each interior boundary loop to the exterior boundary loop
 *  \param prefix the prefix of the output mesh name
 *  the shortest path between the k-th boundary loop and the 0-th boundary loop (exterior one) is 
 *  labeled as sharp edges on the outupt mesh "prefix_k.m".
 */
template <typename M>
void CDynamicShortestPath<M>::_cut( const char *  prefix)
{
  //label all vertices as -1
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
  {
	  M::CVertex * pV = *viter;
	pV->idx() = -1;
  }

  //sort boundary loops by their lengths
	std::vector<M::CLoop*>& loops = m_boundary.loops();

 
  for( size_t k = 0; k < loops.size(); k ++ )
  {
    M::CLoop * pL = loops[k];
	std::list<M::CHalfEdge*>  & pH = pL->halfedges();
    
	for( std::list<M::CHalfEdge*>::iterator hiter = pH.begin(); hiter != pH.end(); hiter ++ )
    {
		M::CHalfEdge * he = *hiter;
		M::CVertex * pV = m_pMesh->halfedgeVertex( he );
        pV->idx() = k;
    }
 }

  //trace the shortest path from the 0-th boundary loop to the k-th boundary loop, output the cut
  for( size_t k = 1; k < loops.size(); k ++ )
  {
      _trace(loops[0], loops[k] );
	//write the sharp edge to the trait string
	 //m_pTrait->write();

	 std::string line;
	 std::stringstream iss(line);
	 iss << prefix << "_" << k-1 << ".cut.vef" ;
	 m_pMesh->write_vef( iss.str().c_str() );
  }

  //reset the edge string
  for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); ++ eiter )
  {
	  M::CEdge * e = *eiter;
	  e->sharp() = false;
  }

  //label sharp edges
  for( std::list<M::CEdge*>::iterator eiter = m_cuts.begin(); eiter != m_cuts.end(); ++ eiter )
  {
	  M::CEdge * e = *eiter;
	  e->sharp() = true;
  }

  //the mesh with all shortest path labeled is outptu to "prefix_cut.m" for computing the fundamental domain
	 std::string line;
	 std::stringstream iss(line);
	 iss << prefix << ".cut.vef" ;
	 m_pMesh->write_vef( iss.str().c_str() );

}

/*!	Compute the shortest path from the source boundary loop to the target boundary loop
 * 
 *  \param source the source boundary loop
 *  \param target the target boundary loop
 */
template<typename M> 
void CDynamicShortestPath<M>::_trace(typename M::CLoop * source , typename M::CLoop * target )
{
	//set all vertex touched to be false
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
    {
		M::CVertex * v = *viter;
		v->touched() = false;
    }

	//set all edge strings to be empty, no edge is sharp
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); ++ eiter )
    {
		M::CEdge * e = *eiter;
		e->sharp() = false;
    }

	std::queue<M::CVertex*> vqueue;
	//enqueue all the target boundary loop vertices	
	for( std::list<M::CHalfEdge * > :: iterator hiter = target->halfedges().begin() ; hiter != target->halfedges().end(); hiter++ )
    {
		M::CHalfEdge * he = *hiter;
		M::CVertex     * pv = m_pMesh->halfedgeVertex( he );
        pv->touched( ) = true;
        vqueue.push( pv );
    }
    
	M::CVertex * destiny = NULL;
	//breadth first search to reach the source loop
	//each vertex has a paraent field indicating the parent vertex
	//in the breadth first searching tree
    while( ! vqueue.empty() )
    {
		M::CVertex * head = vqueue.front();
        vqueue.pop();
        
		for( M::VertexEdgeIterator veiter( head ); !veiter.end(); ++ veiter )
        {
			M::CEdge   * e = *veiter;
			M::CVertex * v1 = m_pMesh->edgeVertex1( e );
			M::CVertex * v2 = m_pMesh->edgeVertex2( e );

			M::CVertex * w = (v1 != head)? v1:v2;

            if( w->idx() == 0 )
            {
                destiny = w;
                w->parent() = head;
				w->bridge() = e;
                break;
            }
           if( w->boundary() ) continue;
            if( w->touched() ) continue;
             w->touched( ) = true;
            vqueue.push( w );
            w->parent() = head;
			w->bridge() = e;
        }
/*
		for( M::VertexEdgeIterator veiter( head ); !veiter.end(); ++ veiter )
        {
			M::CEdge   * e = *veiter;
			M::CVertex * v1 = m_pMesh->edgeVertex1( e );
			M::CVertex * v2 = m_pMesh->edgeVertex2( e );

			M::CVertex * w = (v1 != head)? v1:v2;

            if( w->idx() == 0 )
            {
                destiny = w;
                w->parent() = head;
				w->bridge() = e;
                break;
            }
           if( w->boundary() ) continue;
            if( w->touched() ) continue;
             w->touched( ) = true;
            vqueue.push( w );
            w->parent() = head;
			w->bridge() = e;
        }
*/
        if( destiny != NULL ) break;
    }
	
	//trace back from the source to the target
	//label sharp edges
    assert( destiny != NULL );
	M::CVertex * pv = destiny;
 
    while( pv->parent( ) != NULL )
    {
		M::CVertex * pw = pv->parent( );
		//M::CEdge * e  = m_pMesh->vertexEdge( pv, pw );
		M::CEdge * e  = pv->bridge();
		e->sharp() = true;
		//e->string() = "sharp";
		m_cuts.push_back( e );
		pv = pw;
    }
};


}
#endif