/*!
*      \file DomainCohomology.h
*      \brief Algorithm for Computing Cohomology Group Basis for Mulply Connected Domains
*	   \author David Gu
*      \date Document 10/11/2010
*
*		Computing Cohomology Group Basis for Mulply Connected Domains
*
*/

#ifndef _DOMAIN_COHOMOLOGY_H_
#define _DOMAIN_COHOMOLOGY_H_

#include "CohomologyMesh.h"

namespace MeshLib
{

namespace Topology
{
	/*! \brief CDomainCohomology class
	*  
	*  Algorithm for computing Cohomology Group Basis for Mulply Connected Domains
	*/
  template<typename M>
  class CDomainCohomology
  {
  public:
	  /*! CDomainCohomology constructor
	  *
	  *	\param pMesh input closed mesh
	  * \param pWMesh sliced open mesh
	  */
    CDomainCohomology( M * pMesh, M * pWMesh );
	/*! CHarmonicClosedForm destructor
	*/
    ~CDomainCohomology();
	
	/*!	Compute Cohomology Group Basis for Mulply Connected Domains
	*/
    void calculate_closed_form();

  protected:
	  /*! Input closed mesh */
    M * m_pMesh;
	/*! inpute mesh, sliced open along the shortest paths connecting two boundary loops */
    M * m_pWMesh;

	/*! Compute the closed 1-form*/
	void _closed_1_form();
	/*! Integrate the closed form on the mesh */
	void _integrate();
	
	//label the boundary loop Id on the closed mesh
	/*! label the boundary loop Id on the closed mesh */
	void _label_closed_mesh_boundary_loop_id();
	//locate four corner vertices
	/*! locate four corner vertices 
	*
	*	\param cut_verts store the four corner vertices in the vector cut_verts
	*/
	void _locate_four_corners(  std::vector<typename M::CVertex*> & cut_verts );
	
  };

//CDomainCohomology constructor
//pMesh is the input closed mesh
//pWMesh is the open mesh, obtained by slicing the closed mesh along speical paths
template<typename M>
CDomainCohomology<M>::CDomainCohomology( M * pMesh, M * pWMesh )
{
	m_pMesh  = pMesh; 
	m_pWMesh = pWMesh;
};

//CDomainCohomology destructor
template<typename M>
CDomainCohomology<M>::~CDomainCohomology()
{
};


//Compute the cohomology group basis for multiply connected domains
template<typename M>
void CDomainCohomology<M>::calculate_closed_form()
{
	//compute the closed 1-forms
    _closed_1_form();
	//integrate the harmnoic 1-forms on the mesh
    _integrate();
}





//integrate the closed 1-form on the mesh 
template<typename M>
void CDomainCohomology<M>::_integrate()
{
	for( M::MeshEdgeIterator eiter( m_pWMesh ); !eiter.end(); eiter ++ )
  {
	  M::CEdge * we = *eiter;

	  M::CVertex * w1 = m_pWMesh->edgeVertex1( we );
	  M::CVertex * w2 = m_pWMesh->edgeVertex2( we );

      int id1 = w1->father();
      int id2 = w2->father();

	  M::CVertex * v1 = m_pMesh->idVertex( id1 );
	  M::CVertex * v2 = m_pMesh->idVertex( id2 );

	  M::CEdge * e = m_pMesh->vertexEdge( v1, v2 );

      if( m_pMesh->edgeVertex1( e ) == v1 )
      {
          we->du()  = e->du();
      }
      else
      {  
        we->du() = - e->du();
      }
  }

	M::CVertex * head = NULL;
  for( M::MeshVertexIterator viter( m_pWMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
      v->touched() = false;
      head = v;
  }

  std::queue<M::CVertex*> vqueue;
  head->touched() = true;
  vqueue.push( head );
  head->uv() =  CPoint2(0,0.53);

  while( !vqueue.empty() )
  {
	  M::CVertex * head = vqueue.front();
    vqueue.pop();

	for( M::VertexEdgeIterator veiter(  head ); !veiter.end(); veiter ++ )
    {
		M::CEdge * e = *veiter;
      double du = e->du();

	  M::CVertex * v1 = m_pWMesh->edgeVertex1( e );
	  M::CVertex * v2 = m_pWMesh->edgeVertex2( e );

      if( v1 == head )
      {
        if( v2->touched() ) continue;
        v2->touched( ) = true;
        vqueue.push( v2 );
        v2->uv() = v1->uv() + CPoint2( du,0 );
      }
      else
      {
        assert( v2 == head );
        if( v1->touched() ) continue;
        v1->touched() = true;
        vqueue.push( v1 );
        v1->uv() = v2->uv() - CPoint2( du,0 );
      }
    }
  }
 /* 
  double umin = 1e+10;
  double umax = -1e+10;

  for( CHCFMesh::MeshVertexIterator viter( m_pWMesh );  !viter.end(); viter ++ )
  {
      CHCFVertex* v = *viter;
  	  CPoint2 p = v->uv()
      umax = (umax > p[0] )? umax: p[0];
      umin = (umin < p[0] )? umin: p[0];
  }

 for( MeshVertexIterator viter( m_pWMesh );  !viter.end(); viter ++ )
  {
      CVertex* v = *viter;
  	  CPoint2 p = v_harmonic_closed_form_uv( v );
      p[0] = ( p[0] - umin )/(umax-umin);
      v_harmonic_closed_form_uv( v ) = p;
  }
*/

}


//Compute the closed 1-form
template<typename M>
void CDomainCohomology<M>::_closed_1_form()
{
  for( M::MeshVertexIterator viter( m_pWMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * pV = *viter;
      pV->u() = 0;
  }
 
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * pV = *viter;
      pV->valence() = 0;
  }

  for( M::MeshVertexIterator viter( m_pWMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * pV = *viter;
      int fid = pV->father();
	  M::CVertex * w = m_pMesh->idVertex( fid );
	  w->valence() ++;
  }
  
  
  _label_closed_mesh_boundary_loop_id();

  std::vector<M::CVertex*>  corners;
  _locate_four_corners(  corners );

  //the vertices between corner[1] and corner[2];
  //  2---- 1
  //  |     |
  //  |     |
  //  3---- 0

  std::vector<M::CVertex*> side;
	
    M::CVertex   * pV = corners[1];
	side.push_back( pV );
	M::CVertex   * pN = pV;

	do{
		M::CHalfEdge * pH = m_pWMesh->vertexMostClwOutHalfEdge( pN );
		pN = m_pWMesh->halfedgeVertex( pH );
		side.push_back( pN );
	}while( pN != corners[2] );


 for( int k = 0; k < (int)side.size(); k ++ )
 {
	 M::CVertex * v = side[k];
	 v->u() = 1.0;
 }

 for( M::MeshEdgeIterator eiter( m_pWMesh ); !eiter.end(); eiter ++ )
  {
	  M::CEdge   * we = *eiter;
	  M::CVertex * w1 = m_pWMesh->edgeVertex1( we );
	  M::CVertex * w2 = m_pWMesh->edgeVertex2( we );
	  we->du() = w2->u() - w1->u();

    int id1 =  w1->father();
    int id2 =  w2->father();

	M::CVertex * v1 = m_pMesh->idVertex( id1 );
	M::CVertex * v2 = m_pMesh->idVertex( id2 );

	M::CEdge * e = m_pMesh->vertexEdge( v1, v2 );

    if( m_pMesh->edgeVertex1( e ) == v1 )
    {
        e->du() = we->du();
    }
    else
    {
		assert( m_pMesh->edgeVertex1( e ) == v2 );
        e->du() = -we->du();
    }
  }
}



  //locate 4 corner vertices
  //the vertices between corner[1] and corner[2];
  //  2---- 1
  //  |     |
  //  |     |
  //  3---- 0
template<typename M>
void CDomainCohomology<M>::_locate_four_corners(  std::vector<typename M::CVertex*> & corners )
{

	M::CBoundary boundary( m_pWMesh );
  //the longest one	
	M::CLoop * pL = boundary.loops().front();

	std::list<M::CHalfEdge*> pHs = pL->halfedges();
 
	for( std::list<M::CHalfEdge*>::iterator hiter = pHs.begin(); hiter != pHs.end(); hiter ++)
  {
	  M::CHalfEdge * pH = *hiter;
	  M::CVertex * pV = m_pWMesh->halfedgeVertex( pH );

	  int fid = pV->father();
	  M::CVertex * w = m_pMesh->idVertex( fid );

	  if( !w->boundary()  )		continue;
	  if( w->valence() != 2 )   continue;
	  corners.push_back( pV );
  }

  //cyclicly ajust the order of the vertices in the cut_verts buffer
  //such that the father of the first 2 corners are on the external boundary loop

  assert( corners.size() == 4 );

  std::vector<M::CVertex*> buffer;

  for( int i = 0; i < 4; i ++ )
  {
	  int fid = corners[i]->father();
	  M::CVertex * w = m_pMesh->idVertex( fid );

	  printf("vertex id %d, father id %d, father loop id %d\n", corners[i]->id(), corners[i]->father(),  w->idx() );
  }

  int min_loop_id = 10000;

  for( int i = 0; i < 4; i ++ )
  {
	  M::CVertex * pV = corners[i];
	  M::CVertex * pW  = m_pMesh->idVertex( pV->father() );

	  min_loop_id = ( min_loop_id < pW->idx() )?min_loop_id: pW->idx();
  }

  for( int i = 0; i < 4; i ++ )
  {
	  M::CVertex * pV = corners[(i+0)%4];
	  M::CVertex * pN = corners[(i+1)%4];

	  M::CVertex * pW  = m_pMesh->idVertex( pV->father() );
	  M::CVertex * pWN = m_pMesh->idVertex( pN->father() );

	  if( pW->idx()  == min_loop_id && pWN->idx() == min_loop_id )
	  {
		for( int j = 0; j < 4; j ++ )
		{
			buffer.push_back( corners[(i+j)%4] );
		}
	  }
  }

  corners.clear();

  for( int i = 0; i < 4; i ++ )
  {
	  corners.push_back( buffer[i] );
  }

  for( int i = 0; i < 4; i ++ )
  {
	  int fid = corners[i]->father();
	  M::CVertex * w = m_pMesh->idVertex( fid );

	  printf("vertex id %d, father id %d, father loop id %d\n", corners[i]->id(), corners[i]->father(),  w->idx() );
  }

};

//label closed mesh boundary ids
template<typename M>
void CDomainCohomology<M>::_label_closed_mesh_boundary_loop_id()
{
  //sort boundary loops by their lengths
  
	M::CBoundary boundary( m_pMesh );

	std::vector<M::CLoop*>& loops = boundary.loops();

	
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

};



} //namespace Topology
} //namespace MeshLib

#endif _DOMAIN_COHOMOLOGY_H_