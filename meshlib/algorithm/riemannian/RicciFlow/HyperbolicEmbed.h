/*! \file HyperbolicEmbed.h
 *  \brief Isometrically embed a mesh with hyperbolic metric onto the Poincare disk
 *  \author David Gu
 *  \date  documented on 10/17/2010
 *
 *  Embed a mesh with hyperbolic metric on the Poincare disk
 */
#ifndef  _HYPERBOLIC_EMBED_H_
#define  _HYPERBOLIC_EMBED_H_

#include <vector>
#include <queue>

#include "RicciFlowMesh.h"
#include "BaseEmbed.h"
#include "Geometry/HyperbolicCircle.h"

namespace MeshLib
{

namespace RicciFlow
{

/*! \brief CHyperbolicEmbed class
 *
 *   embed a mesh with flat metric onto the Poincare disk
 */
template<typename M>
class CHyperbolicEmbed : public CBaseEmbed<M>
  {
  public: 
   /*! \brief CHyperbolicEmbed constructor */
    CHyperbolicEmbed( M * pMesh );
   /*! \brief CHyperbolicEmbed destructor */
    ~CHyperbolicEmbed(){};
	/*! _embed the mesh */
    void _embed();
  protected:
	/*! embed the first face 
	 * \param head the first face
	 */
	  void _embed_first_face( typename M::CFace * head );
	/*! initialization */
	void _initialize();
	/*! embed one face
	 *  \param face the face to be embed */
	void _embed_face( typename M::CFace * face  );

  };


//constructor
template<typename M>
CHyperbolicEmbed<M>::CHyperbolicEmbed( M * pMesh ): CBaseEmbed( pMesh )
{
};

//embed the first face
template<typename M>
void CHyperbolicEmbed<M>::_embed_first_face( typename M::CFace * head )
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

	av[0]->huv() = CPoint2(0,0 );
	
	M::CEdge * e = m_pMesh->halfedgeEdge( he[0] );
	double r = e->length();
	av[1]->huv() = CPoint2( (exp(r)-1.0)/(exp(r)+1.0), 0 );
	
	r = m_pMesh->halfedgeEdge( he[2] )->length();
	r = (exp(r)-1.0)/(exp(r)+1.0);

	double l[3];
	for( int i = 0; i < 3; i ++ )
	{
		M::CEdge * e = m_pMesh->halfedgeEdge( he[i] );
	  l[i] = e->length();
	}
	
	double phi = acos( (cosh(l[2]) * cosh(l[0])-cosh(l[1]))/(sinh(l[2])*sinh(l[0])) );

	double cs = cos( phi );
	double sn = sin( phi );

	av[2]->huv() = CPoint2(cs *r, sn * r );

	
	//set touched field to be true

	for(int i = 0; i < 3; i ++ )
	{
		m_pMesh->halfedgeTarget( he[i] )->touched() = true;
	}

};




template<typename M>
void CHyperbolicEmbed<M>::_initialize()
{
	for( M::MeshVertexIterator viter(  m_pMesh ); !viter.end(); viter ++)
  {
	  M::CVertex * v = *viter;
    v->huv() = CPoint2( 0, 0);
    v->touched() = false;
  }


  for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
  {
	  M::CFace * f = *fiter;
    f->touched() = false;
  }

  M::CFace * root_face = NULL;
  for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
  {
	  M::CFace * f = *fiter;
	int sharp_count = 0;
	for( M::FaceHalfedgeIterator fhiter( f ); !fhiter.end(); fhiter ++ )
	{
		M::CHalfEdge * h = *fhiter;
		M::CEdge * e = m_pMesh->halfedgeEdge( h );
		if( e->sharp() ) sharp_count ++;
	}
	if( sharp_count == 3 ) 
	{
		root_face = f;
		break;
	}
  }

  if( root_face == NULL )
  {
	  printf("Error: You need to use sharp edges to specify the root_face.\n");
	  exit(-1);
  }
  	
	_embed_first_face( root_face );
	root_face->touched() = true;

	M::CHalfEdge * he = m_pMesh->faceMostCcwHalfEdge( root_face );

  do{
	  M::CHalfEdge  * sh = m_pMesh->halfedgeSym( he );
    if( sh != NULL )
    {
		M::CFace * df = m_pMesh->halfedgeFace( sh );
      m_queue.push( df );
      df->touched() = true;
    }
    he = m_pMesh->faceNextCcwHalfEdge( he );
  }while ( he != m_pMesh->faceMostCcwHalfEdge( root_face ) );

 };


template<typename M>
void CHyperbolicEmbed<M>::_embed()
{

	_initialize();

	int count = 0;

	while( !m_queue.empty() )
	{
		M::CFace *  head = m_queue.front();
		m_queue.pop();
		assert( head->touched() );
		_embed_face( head );

		M::CHalfEdge * he = m_pMesh->faceMostCcwHalfEdge( head );
		do{
			M::CHalfEdge * sh = m_pMesh->halfedgeSym( he );
		  if( sh != NULL )
		  {  
			  M::CFace * df = m_pMesh->halfedgeFace( sh );
			if( !df->touched() )
			{
				m_queue.push( df );
				df->touched() = true;
			}
		  }
		  he = m_pMesh->halfedgeNext( he );
		}while( he != m_pMesh->faceMostCcwHalfEdge( head ) );
	}
};


template<typename M>
void CHyperbolicEmbed<M>::_embed_face( typename M::CFace * head )
{

	std::vector<M::CVertex*> av;
  for( M::FaceVertexIterator fviter( head ); !fviter.end(); fviter ++ )
  {
	  M::CVertex * pV = *fviter;
	  av.push_back( pV );
  }

  M::CVertex * A = NULL;
  M::CVertex * B = NULL;
  M::CVertex * C = NULL;

	for( int i = 0; i < 3; i ++ )
	{
	  if( av[i]->touched() ) continue;
	  C = av[(i+0)%3];
	  A = av[(i+1)%3];
	  B = av[(i+2)%3];
	  break;
	}
		
	if( C == NULL ) return;

	assert( A->touched() && B->touched() );

	//radius of the first circle
	M::CEdge * e = m_pMesh->vertexEdge(A,C);
	double radius_1 = e->length();

	e = m_pMesh->vertexEdge(B,C);
	//radius of the second circle
	double radius_2 = e->length();

	//center of the first circle
	CPoint2 c1 = A->huv();

	//center of the second circle
	CPoint2 c2 = B->huv();

	CHyperbolicCircle C1( c1, radius_1 );
	CHyperbolicCircle C2( c2, radius_2 );

	CPoint2 i1;
	CPoint2 i2;



	int result = C1.intersect( C2, i1, i2 );


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
}


typedef CHyperbolicEmbed<CRFMesh> CHRFEmbed;	

} //namespace RicciFlow
} //namespace MeshLib

#endif 