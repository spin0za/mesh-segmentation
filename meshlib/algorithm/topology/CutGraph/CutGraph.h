/*!
*      \file CutGraph.h
*      \brief Algorithm for Cut Graph
*	   \author David Gu
*      \date Document 12/24/2010
*
*		Computing cut graph
*
*/

/********************************************************************************************************************************
*
*      Cut Graph Class
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       Compute Cut Graph
* 
*       David Gu December 24, 2010,  gu@cs.stonybrook.edu
*
*
*      Input:
*         
*           Closed Surface
*
*      Output:
*
*           Cut graph on the surface labelled as sharp edges
*
*********************************************************************************************************************************/

/*---------------------------------------------------------------------------------------------------------------------------------
#include <math.h>
#include "mesh/mesh.h"
#include "CutGraph/CutGraph.h"

using namespace MeshLib;

int main( int argc, char * argv[] )
{

  if( strcmp( argv[1], "-cut_graph" ) == 0 )
  {
	CutGraphMesh mesh;
	mesh.read_m( argv[2] ); //original mesh

	CCutGraph cut_graph( &mesh );
	cut_graph.propagate();

	mesh.write_m( argv[3] );
	return 0;
  }

}
----------------------------------------------------------------------------------------------------------------------------------*/

#ifndef _CUT_GRAPH_H_
#define _CUT_GRAPH_H_

#include  <math.h>
#include <queue>
#include "CutGraphMesh.h"

namespace MeshLib
{

namespace Topology
{
	/*! \brief CCutGraph class
	*  
	*  Algorithm for computing cut graph
	*/
  template<typename M>
  class CCutGraph
  {
  public:
	  /*! CCutGraph constructor
	  *
	  *	\param pMesh input closed mesh
	  */

    CCutGraph( M * pMesh  );
	/*!
	 * CCutGraph destructor
	 */
    ~CCutGraph();
	/*!
	 *	Compute the spanning tree of the dual mesh
	 *  the edges whose duals are not on the tree form the cut locus
	 *  label the cut locus as the sharp edges
	 */
	void _cut_locus();
	/*!
	 *	Compute the spanning tree of the dual mesh
	 *  the edges whose duals are not on the tree form the cut locus
	 *  label the cut locus as the sharp edges
	 */
	void _cut_locus( int head_face_id );

  protected:
    /*!
	 *	Input closed mesh
	 */ 
    CutGraphMesh * m_pMesh;

	/*!
	 *	Compute the spanning tree of the dual mesh
	 *  the edges whose duals are not on the tree are labelled as the sharp edges
	 */
	void _dual_spanning_tree();

	/*!
	 *	Compute the spanning tree of the dual mesh
	 *  the edges whose duals are not on the tree are labelled as the sharp edges
	 */
	void _dual_spanning_tree( int head_face_id );


	/*!
	 *	Prune the cut locus
	 */
    void _prune();

  };

template<typename M>
CCutGraph<M>::CCutGraph( M * pMesh )
{
	m_pMesh = pMesh;
}

template<typename M>
CCutGraph<M>::~CCutGraph()
{
}

template<typename M>
void CCutGraph<M>::_dual_spanning_tree()
{
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		e->sharp()= false;
	}

	M::CFace * head = NULL;
	for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		M::CFace * f = *fiter;
		f->touched() = false;
		head = f;
	}

  assert( head != NULL );

  head->touched() = true;
  
  std::queue<M::CFace*>  fqueue;
  fqueue.push( head );

  while( !fqueue.empty() )
  {
	  M::CFace * head = fqueue.front();
    fqueue.pop();

	for( M::FaceHalfedgeIterator fhiter( head ); !fhiter.end(); fhiter ++ )
	{
		M::CHalfEdge * he = *fhiter;
		M::CHalfEdge * sh = m_pMesh->halfedgeSym( he );

		  if( sh != NULL )
		  {
			M::CFace * sf = m_pMesh->halfedgeFace( sh );
			if( !sf->touched() )
			{
				fqueue.push( sf );
				sf->touched() = true;
				M::CEdge * pE = m_pMesh->halfedgeEdge( he );
				pE->sharp() = true;
			}
		  }
    }

  }

  for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
  {
	  M::CEdge * e = *eiter;
		e->sharp() = !e->sharp();
  }

}

template<typename M>
void CCutGraph<M>::_cut_locus()
{
	_dual_spanning_tree();
	_prune();
}

template<typename M>
void CCutGraph<M>::_prune()
{
	//compute the spanning tree of the cut locus
	std::set<M::CVertex*> nodes;

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pv = *viter;
		pv->valence() = 0;
		pv->touched() = false;
	}

	M::CVertex * head = NULL;
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		e->mark() = 0;

		if( !e->sharp() ) continue;
		e->mark() = 1;
		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		head = v1;
	}

	M::CVertex * root = head;

	std::queue<M::CVertex*> vqueue;
	head->touched() = true;
	vqueue.push( head );
	head->father() = NULL;

	while( !vqueue.empty() )
	{
		head = vqueue.front();
		vqueue.pop();

		for( M::VertexEdgeIterator veiter( head ); !veiter.end(); veiter ++ )
		{
			M::CEdge * e = *veiter;
			if( e->mark() == 1 )
			{
				M::CVertex * v1 = m_pMesh->edgeVertex1( e );
				M::CVertex * v2 = m_pMesh->edgeVertex2( e );
				M::CVertex * tail = (v1 != head )?v1:v2;
				if( tail->touched() ) continue;
				vqueue.push( tail );
				tail->touched() = true;
				tail->father()  = head;
				e->mark() = 2;
			}
		}
	}

    //compute the edges on the cut locus, but not on the spanning tree

	std::list<M::CEdge*> ecuts;

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter; 
		if( e->mark() == 1 )
		{
			ecuts.push_back( e );
		}
		e->sharp() = false;
	}

	//each cut edge form a loop on the cut locus

	for( std::list<M::CEdge*>::iterator eiter = ecuts.begin(); eiter != ecuts.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		e->sharp() = true;
		
		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );

		M::CVertex * pV = v1;
		while( pV->father() != NULL )
		{
			M::CVertex * fv = pV->father();
			M::CEdge * e = m_pMesh->vertexEdge( pV, fv );
			e->sharp() = true;
			pV = fv;
		}
		pV = v2;
		while( pV->father() != NULL )
		{
			M::CVertex * fv = pV->father( );
			M::CEdge * e = m_pMesh->vertexEdge( pV, fv );
			e->sharp() = true;
			pV = fv;
		}
	}

	//the root might be with valence 1, then remove the link attached to the root
	while( true )
	{
		M::CEdge * link = NULL;
		int valence = 0;
		for( M::VertexEdgeIterator veiter( root ); !veiter.end(); ++ veiter )
		{
			M::CEdge * e = *veiter;
			if( !e->sharp() ) continue;
			valence ++;
			link = e;
		}
		if( valence != 1 ) break;
		link->sharp( ) = false;

		M::CVertex * v1 = m_pMesh->edgeVertex1( link );
		M::CVertex * v2 = m_pMesh->edgeVertex2( link );

		root = ( root != v1 )?v1:v2;
	}

}

template<typename M>
void CCutGraph<M>::_cut_locus( int head_face_id )
{
	_dual_spanning_tree( head_face_id );
	_prune();
}

template<typename M>
void CCutGraph<M>::_dual_spanning_tree( int head_face_id )
{
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		e->sharp()= false;
	}

	M::CFace * head = m_pMesh->idFace( head_face_id );

	assert( head != NULL );

  head->touched() = true;
  
  std::queue<M::CFace*>  fqueue;
  fqueue.push( head );

  while( !fqueue.empty() )
  {
	  M::CFace * head = fqueue.front();
    fqueue.pop();

	for( M::FaceHalfedgeIterator fhiter( head ); !fhiter.end(); fhiter ++ )
	{
		M::CHalfEdge * he = *fhiter;
		M::CHalfEdge * sh = m_pMesh->halfedgeSym( he );

		  if( sh != NULL )
		  {
			  M::CFace * sf = m_pMesh->halfedgeFace( sh );
			if( !sf->touched() )
			{
				fqueue.push( sf );
				sf->touched() = true;
				M::CEdge * pE = m_pMesh->halfedgeEdge( he );
				pE->sharp() = true;
			}
		  }
    }

  }

  for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
  {
	  M::CEdge * e = *eiter;
		e->sharp() = !e->sharp();
  }

}

} //namespace Topology
} //MeshLib
#endif _CUT_MESH_H_