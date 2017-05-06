/*!
*      \file Homology.h
*      \brief Algorithm for Homology
*	   \author David Gu
*      \date Document 12/25/2010
*
*		Computing Homology
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
*       Compute Homology Group Basis
* 
*       David Gu December 25, 2010,  gu@cs.stonybrook.edu
*
*
*      Input:
*         
*           Closed Surface
*
*      Output:
*
*           Homology group basis, labelled as sharp edges
*
*********************************************************************************************************************************/

/*---------------------------------------------------------------------------------------------------------------------------------
#include <math.h>
#include "mesh/mesh.h"
#include "CutGraph/Homology.h"

using namespace MeshLib;

int main( int argc, char * argv[] )
{
	if( strcmp( argv[1] , "-homology") == 0 )
	{
		CutGraphMesh spm;
		spm.read_m( argv[2] );

		CHomology homology( & spm );
		homology._calculate_base();
		homology._output( argv[3] );
		return 0;
	}
}
----------------------------------------------------------------------------------------------------------------------------------*/

#ifndef _HOMOLOGY_H_
#define _HOMOLOGY_H_

#include <queue>
#include <vector>

#include "CutGraph.h"
#include "Graph.h"


namespace MeshLib
{

namespace Topology
{
/*! \brief CHomology class
*  
*  Algorithm for computing Homology Group Basis
*/
template<typename M>
class	CHomology
{
public:
	/*!	CHomology Class Constructor
	 *
	 *  \param pMesh input closed mesh
	 */
	CHomology( M * pMesh );
	/*!	CHomology Class destructor
	 */
	~CHomology();
	/*!	Compute the homology group basis
	 */
	void _calculate_base();
	/*!	Output each homology basis to a file, the loop is labelled as sharp edges
	 */
	void _output( const char * prefix );
	
protected:
	/*!
	 *	Input closed mesh
	 */
	M * m_pMesh;
	/*! Cut graph
	 */
	CCutGraph<M>    * m_cut_graph;
	/*! Compute the loop basis of the cut graph
	 */
	CGraph<M>        * m_graph;
	/*!	The list of homology base loop
	 */
	std::list<std::list<typename M::CHalfEdge*>*> m_loops;
};

template<typename M>
CHomology<M>::CHomology( M * pMesh )
{
    m_pMesh = pMesh;

    m_cut_graph = new CCutGraph<M>( m_pMesh );
    m_graph = new CGraph<M>( m_pMesh );
}

template<typename M>
CHomology<M>::~CHomology()
{
  delete m_cut_graph;
  delete m_graph;
    
  for( std::list<std::list<M::CHalfEdge*>*>::iterator liter = m_loops.begin(); liter != m_loops.end(); liter ++ )
    {
		std::list<M::CHalfEdge*>* pL = *liter;
        delete pL;
     }
    m_loops.clear();
}

template<typename M>
void CHomology<M>::_calculate_base()
{
    m_cut_graph->_cut_locus();

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
    {
		M::CEdge * e = *eiter;
	  if( e->boundary()) continue; //we do not consider boundary loops, 3/22/2011
      if( !e->sharp() ) continue;

      m_graph->_insert( e );
    }

    m_graph->_breadth_first_search();
    m_graph->_locate_loops();
    
	std::list<std::list<M::CHalfEdge*>*> & loops = m_graph->loops();

	for( std::list<std::list<M::CHalfEdge*>*>::iterator liter = loops.begin(); liter != loops.end(); liter ++ )
    {
		std::list<M::CHalfEdge*>* pL = *liter;
		std::list<M::CHalfEdge*>* pLoop = new std::list<M::CHalfEdge*>;
        assert( pLoop );
        (*pLoop) = (*pL);
        m_loops.push_back( pLoop );
    }
}
template<typename M>
void CHomology<M>::_output( const char * prefix )
{
  int id = 0;
  for( std::list<std::list<M::CHalfEdge*>*>::iterator liter = m_loops.begin(); liter != m_loops.end(); liter ++ )
  {

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * pE = *eiter;
		pE->sharp() = false;
	}

	std::list<M::CHalfEdge*>* pL = *liter;
    
	for( std::list<M::CHalfEdge*>::iterator hiter = pL->begin(); hiter != pL->end(); hiter ++ )
	{
		M::CHalfEdge * he = *hiter;
		M::CEdge * e = m_pMesh->halfedgeEdge(he);
		e->sharp() = true;
	}

	std::stringstream iss;
	iss << prefix << "_" << id++ << ".m" ;
	m_pMesh->write_m( iss.str().c_str() );
  }
};


} //namespace Topology
} //namespace MeshLib
#endif  //_HOMOLOGY_H_