/*! \file EuclideanEmbed.h
 *  \brief Isometrically embed a mesh with flat metric onto the plane
 *  \author David Gu
 *  \date  documented on 10/17/2010
 *
 *  Embed a mesh with flat metric on canonical domain
 */
#ifndef  _EUCLIDEAN_EMBED_H_
#define  _EUCLIDEAN_EMBED_H_

#include <vector>
#include <queue>

#include "RicciFlowMesh.h"
#include "BaseEmbed.h"
#include "Geometry/Circle.h"

namespace MeshLib
{

namespace RicciFlow
{

/*! \brief CEmbed class
 *
 *   embed a mesh with flat metric onto the plane
 */
template<typename M>
class CEuclideanEmbed : public CBaseEmbed<M>
  {
  public: 
   /*! \brief CEuclideanEmbed constructor */
    CEuclideanEmbed( M * pMesh );
   /*! \brief CEuclideanEmbed destructor */
    ~CEuclideanEmbed(){};

  protected:
	/*! embed the first face 
	 * \param head the first face
	 */
	  void _embed_first_face( typename M::CFace * head );
	/*!
	 *	embed one face
	 */
	  void _embed_face( typename M::CFace * head );

  };


//constructor
template<typename M>
CEuclideanEmbed<M>::CEuclideanEmbed( M * pMesh ):CBaseEmbed( pMesh )
{
};

//embed the first face
template<typename M>
void CEuclideanEmbed<M>::_embed_first_face( typename M::CFace * head )
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

template<typename M>
void CEuclideanEmbed<M>::_embed_face( typename M::CFace * head )
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

	//radius of the first circle
	double r1 = m_pMesh->vertexEdge(A,C)->length();
	//radius of the second circle
	double r2 = m_pMesh->vertexEdge(B,C)->length();
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
}


typedef CEuclideanEmbed<CRFMesh> CRFEmbed;	

}	//namespace RicciFlow
}	//namespace MeshLib


#endif 