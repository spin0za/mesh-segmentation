/*!
*      \file ShortestPath.cpp
*      \brief Implement CShortestPath class
*	   \author David Gu
*      \date Document 10/11/2010
*
*		Compute the shortest path between two boundary loops, and labeled as sharp edges
*
*/
#include "ShortestPath.h"

using namespace MeshLib;

/*! CShortestPath constructor
 *
 *  \param pMesh the input mesh
 */

CShortestPath::CShortestPath( CSPMesh * pMesh ): m_pMesh( pMesh ), m_boundary( m_pMesh )
{
}

/*! CShortestPath destructor
 *
 */


CShortestPath::~CShortestPath()
{
}

/*!	Compute the shortest paths between each interior boundary loop to the exterior boundary loop
 *  \param prefix the prefix of the output mesh name
 *  the shortest path between the k-th boundary loop and the 0-th boundary loop (exterior one) is 
 *  labeled as sharp edges on the outupt mesh "prefix_k.m".
 */
void CShortestPath::_cut( const char *  prefix)
{
  //label all vertices as -1
	for( CSPMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
  {
    CSPVertex * pV = *viter;
	pV->idx() = -1;
  }

  //sort boundary loops by their lengths
	std::vector<CSPMesh::CLoop*>& loops = m_boundary.loops();

 
  for( size_t k = 0; k < loops.size(); k ++ )
  {
    CSPMesh::CLoop * pL = loops[k];
    std::list<CHalfEdge*>  & pH = pL->halfedges();
    
    for( std::list<CHalfEdge*>::iterator hiter = pH.begin(); hiter != pH.end(); hiter ++ )
    {
        CHalfEdge * he = *hiter;
        CSPVertex * pV = m_pMesh->halfedgeVertex( he );
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
	 iss << prefix << "_" << k-1 << ".cut.m" ;
	 m_pMesh->write_m( iss.str().c_str() );
  }

  //reset the edge string
  for( CSPMesh::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); ++ eiter )
  {
	  CSPEdge * e = *eiter;
	  e->string() = "";
	  e->sharp() = false;
  }

  //label sharp edges
  for( std::list<CSPEdge*>::iterator eiter = m_cuts.begin(); eiter != m_cuts.end(); ++ eiter )
  {
	  CSPEdge * e = *eiter;
	  e->sharp() = true;
	  e->string() = "sharp";
  }

  //the mesh with all shortest path labeled is outptu to "prefix_cut.m" for computing the fundamental domain
	 std::string line;
	 std::stringstream iss(line);
	 iss << prefix << ".cut.m" ;
	 m_pMesh->write_m( iss.str().c_str() );

}

/*!	Compute the shortest path from the source boundary loop to the target boundary loop
 * 
 *  \param source the source boundary loop
 *  \param target the target boundary loop
 */

void CShortestPath::_trace(CSPMesh::CLoop * source , CSPMesh::CLoop * target )
{
	//set all vertex touched to be false
	for( CSPMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
    {
      CSPVertex * v = *viter;
      v->touched() = false;
    }

	//set all edge strings to be empty, no edge is sharp
	for( CSPMesh::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); ++ eiter )
    {
      CSPEdge * e = *eiter;
      e->sharp() = false;
      e->string() = "";
    }

    std::queue<CSPVertex*> vqueue;
	//enqueue all the target boundary loop vertices	
    for( std::list<CHalfEdge * > :: iterator hiter = target->halfedges().begin() ; hiter != target->halfedges().end(); hiter++ )
    {
        CHalfEdge * he = *hiter;
        CSPVertex     * pv = m_pMesh->halfedgeVertex( he );
        pv->touched( ) = true;
        vqueue.push( pv );
    }
    
    CSPVertex * destiny = NULL;
	//breadth first search to reach the source loop
	//each vertex has a paraent field indicating the parent vertex
	//in the breadth first searching tree
    while( ! vqueue.empty() )
    {
		CSPVertex * head = vqueue.front();
        vqueue.pop();
        
		for( CSPMesh::VertexVertexIterator vviter( head ); !vviter.end(); ++ vviter )
        {
            CSPVertex * w = *vviter;
            if( w->idx() == 0 )
            {
                destiny = w;
                w->parent() = head;
                break;
            }
           if( w->boundary() ) continue;
            if( w->touched() ) continue;
             w->touched( ) = true;
            vqueue.push( w );
            w->parent() = head;
        }

        if( destiny != NULL ) break;
    }
	
	//trace back from the source to the target
	//label sharp edges
    assert( destiny != NULL );
    CSPVertex * pv = destiny;
 
    while( pv->parent( ) != NULL )
    {
      CSPVertex * pw = pv->parent( );
      CSPEdge * e  = m_pMesh->vertexEdge( pv, pw );
      e->sharp() = true;
	  e->string() = "sharp";
	  m_cuts.push_back( e );
      pv = pw;
    }
}

