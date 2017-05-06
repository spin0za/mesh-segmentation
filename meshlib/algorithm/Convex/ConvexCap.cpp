/*!
*      \file ConvexCap.cpp
*      \brief Implement CConvexCap class
*	   \author David Gu
*      \date Document 03/22/2013
*
*		isometric embedding of convex cap
*
*/


#include "ConvexCap.h"
#include "Structure/Structure.h"
#include "Mesh/iterators.h"

using namespace MeshLib;
using std::vector;

//Constructor
//Count the number of interior vertices
//Count the number of boundary vertices
//Compute edge weight

/*!	CHarmonicMapper constructor 
*	Count the number of interior vertices, boundary vertices and the edge weight
*
*/
CConvexCap::CConvexCap( CCCMesh* pMesh ): m_pMesh( pMesh ), m_boundary( m_pMesh )
{
	for( CCCMesh::MeshEdgeIterator eiter( m_pMesh); !eiter.end(); eiter ++ )
	{
		CConvexCapEdge * pE = *eiter;
		pE->source_l() = m_pMesh->edgeLength( pE );
	}

	CPoint max(-1e+10,-1e+10,0);
	CPoint min( 1e+10, 1e+10,0);

	for( CCCMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		CConvexCapVertex * pV = *viter;
		CPoint p = pV->point();

		for( int i = 0; i < 2; i ++ )
		{
			max[i] =  (max[i]>p[i])? max[i]:p[i];
			min[i] =  (min[i]<p[i])? min[i]:p[i];
		}
	}

	for( CCCMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		CConvexCapVertex * pV = *viter;
		pV->point() -= (max+min)/2.0;
	}

	double s = (max[0]-min[0]>max[1]-min[1])?max[0]-min[0]:max[1]-min[1];
	s/= 2.0;

	for( CCCMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		CConvexCapVertex * pV = *viter;
		pV->point() /= s;
	}

	for( CCCMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		CConvexCapVertex * pV = *viter;
		if( m_pMesh->isBoundary( pV ) )
		{
			pV->point()[2] = 0;
			continue;
		}
		CPoint p = pV->point();
		pV->point()[2] = 0.05 * (1.0 - p[0] * p[0] - p[1] * p[1]);
	}

	m_pMesh->write_m("test.m");

	for( CCCMesh::MeshEdgeIterator eiter( m_pMesh); !eiter.end(); eiter ++ )
	{
		CConvexCapEdge * pE = *eiter;
		pE->target_l() = m_pMesh->edgeLength( pE );
	}

	for( CCCMesh::MeshEdgeIterator eiter( m_pMesh); !eiter.end(); eiter ++ )
	{
		CConvexCapEdge * pE = *eiter;
		//pE->l() = 0.95 * pE->source_l() + 0.05 * pE->target_l();
		pE->l() = pE->target_l();
	}

	for( CCCMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		CConvexCapVertex * pV = *viter;
		if( m_pMesh->isBoundary( pV ) )
		{
			pV->h() = 0;
			continue;
		}
		//pV->h() = 1e-3;
		pV->h()  = 0;
	}

	_create_interior_prism();
	_create_boundary_prism();
}

void CConvexCap::_create_boundary_prism()
{
	
	for(int i = 0; i < 4; i ++ )
	{
		CConvexCapVertex * pv = m_boundary_prism.createVertex(i);
	}

	std::vector<CConvexCapVertex*> vs;

	vs.push_back( m_boundary_prism.idVertex(0) );
	vs.push_back( m_boundary_prism.idVertex(2) );
	vs.push_back( m_boundary_prism.idVertex(1) );

	m_boundary_prism.createFace( vs, 0 );

	vs.clear();
	vs.push_back( m_boundary_prism.idVertex(1) );
	vs.push_back( m_boundary_prism.idVertex(2) );
	vs.push_back( m_boundary_prism.idVertex(3) );

	m_boundary_prism.createFace( vs, 1 );

	vs.clear();
	vs.push_back( m_boundary_prism.idVertex(2) );
	vs.push_back( m_boundary_prism.idVertex(0) );
	vs.push_back( m_boundary_prism.idVertex(3) );

	m_boundary_prism.createFace( vs, 2 );

	vs.clear();
	vs.push_back( m_boundary_prism.idVertex(0) );
	vs.push_back( m_boundary_prism.idVertex(1) );
	vs.push_back( m_boundary_prism.idVertex(3) );

	m_boundary_prism.createFace( vs, 3 );

}


void CConvexCap::_create_interior_prism()
{
	
	for(int i = 0; i < 6; i ++ )
	{
		CConvexCapVertex * pv = m_interior_prism.createVertex(i);
	}

	std::vector<CConvexCapVertex*> vs;

	vs.push_back( m_interior_prism.idVertex(0) );
	vs.push_back( m_interior_prism.idVertex(2) );
	vs.push_back( m_interior_prism.idVertex(1) );

	m_interior_prism.createFace( vs, 0 );

	vs.clear();
	vs.push_back( m_interior_prism.idVertex(3) );
	vs.push_back( m_interior_prism.idVertex(4) );
	vs.push_back( m_interior_prism.idVertex(5) );

	m_interior_prism.createFace( vs, 1 );

	vs.clear();
	vs.push_back( m_interior_prism.idVertex(0) );
	vs.push_back( m_interior_prism.idVertex(1) );
	vs.push_back( m_interior_prism.idVertex(4) );
	vs.push_back( m_interior_prism.idVertex(3) );

	m_interior_prism.createFace( vs, 2 );

	vs.clear();
	vs.push_back( m_interior_prism.idVertex(1) );
	vs.push_back( m_interior_prism.idVertex(2) );
	vs.push_back( m_interior_prism.idVertex(5) );
	vs.push_back( m_interior_prism.idVertex(4) );

	m_interior_prism.createFace( vs, 3 );

	vs.clear();
	vs.push_back( m_interior_prism.idVertex(2) );
	vs.push_back( m_interior_prism.idVertex(0) );
	vs.push_back( m_interior_prism.idVertex(3) );
	vs.push_back( m_interior_prism.idVertex(5) );

	m_interior_prism.createFace( vs, 4 );

}



void CConvexCap::__compute_prism_face_normal( CCCMesh & _prism )
{
	for( CCCMesh::MeshFaceIterator fiter( &_prism ); !fiter.end(); fiter ++ )
	{
		CConvexCapFace * pF = *fiter;
		std::vector<CConvexCapVertex*> vs;
		for( CCCMesh::FaceVertexIterator fviter( pF ); !fviter.end(); fviter ++ )
		{
			CConvexCapVertex * pV = *fviter;
			vs.push_back( pV );
		}
	}
}

/*!
 *	Compute edge dihedral angle of prism
 */
void CConvexCap::__compute_prism_edge_dihedral_angle( CCCMesh & _prism )
{
	for( CCCMesh::MeshEdgeIterator eiter( &_prism ); !eiter.end(); eiter ++ )
	{
		CConvexCapEdge * pE = *eiter;

		CConvexCapHalfEdge * ph0 = _prism.edgeHalfedge(pE,0);
		CConvexCapHalfEdge * ph1 = _prism.edgeHalfedge(pE,1);

		CConvexCapFace * pF0 = _prism.halfedgeFace( ph0 );
		CConvexCapFace * pF1 = _prism.halfedgeFace( ph1 );
		
		CConvexCapVertex * pT = _prism.halfedgeTarget( ph0 );
		CConvexCapVertex * pS = _prism.halfedgeSource( ph0 );

		CPoint d = pT->point() - pS->point();
		d = d/d.norm();

		CPoint c = pF0->normal()^pF1->normal();
		double ang = asin( c * d);
		if( ang < 0 ) ang += 3.1415926535;

		pE->angle()= ang;
	}
}

//Destructor
/*!
 *	CConvexCap destructor
 */
CConvexCap::~CConvexCap()
{
}


void CConvexCap::_embed_interior_prism( CConvexCapFace * pFace, CCCMesh & _prism )
{
	std::vector<CConvexCapHalfEdge*> hes;
	std::vector<CConvexCapVertex*> vs;

	for( CCCMesh::FaceHalfedgeIterator fhiter( pFace ); !fhiter.end(); fhiter ++ )
	{
		CConvexCapHalfEdge * pH = *fhiter;
		hes.push_back( pH );
		CConvexCapVertex * pV = m_pMesh->halfedgeSource( pH ); 
		vs.push_back( pV );
	}

	/* compute the length */
	std::vector<double> ls;

	for( size_t i = 0; i < hes.size(); i ++ )
	{
		CConvexCapHalfEdge * pH = hes[i];
		CConvexCapEdge     * pE = m_pMesh->halfedgeEdge( pH );		
		CConvexCapVertex * pT = m_pMesh->halfedgeTarget( pH );
		CConvexCapVertex * pS = m_pMesh->halfedgeSource( pH );
		
		double  l = pE->l();
		double hs = pS->h();
		double ht = pT->h();

		l = sqrt( l*l - (hs-ht) * (hs - ht) );
		ls.push_back( l );
	}

	std::vector<double> omega;
	for( int i = 0; i < 3; i ++ )
	{
		double d = ls[(i+1)%3]*ls[(i+1)%3] + ls[(i+2)%3] * ls[(i+2)%3] - ls[(i+0)%3] * ls[(i+0)%3];
		d = d/( 2 * ls[(i+1)%3] * ls[(i+2)%3] );
		omega.push_back( acos( d ) );
	}

	


	for( int i = 0; i < 3; i ++ )
	{
		hes[i]->omega() = omega[(i+2)%3];
	}

	CConvexCapVertex * pV = _prism.idVertex(0);
	pV->point() = CPoint(0,0,0);
	pV = _prism.idVertex(1);
	pV->point() = CPoint(ls[0],0,0);
	pV = _prism.idVertex(2);
	pV->point() = CPoint(ls[2]*cos(omega[1]),ls[2]*sin(omega[1]),0);

	for( int i = 0; i < 3; i ++ )
	{
		CConvexCapVertex * pV = _prism.idVertex(i+0);
		CConvexCapVertex * pW = _prism.idVertex(i+3);
		pW->point() = pV->point() + CPoint(0,0, vs[i]->h() );
	}
}


void CConvexCap::_embed_boundary_prism( CConvexCapFace * pFace, CCCMesh & _prism )
{
	std::vector<CConvexCapHalfEdge*> hs;

	for( CCCMesh::FaceHalfedgeIterator fhiter( pFace ); !fhiter.end(); fhiter ++ )
	{
		CConvexCapHalfEdge * pH = *fhiter;
		hs.push_back( pH );
	}
	//set the first halfedge to be the boundary halfedge
	int start = 0;
	for( size_t i = 0; i < 3; i++ )
	{
		CConvexCapHalfEdge * pH = hs[i];
		CConvexCapEdge     * pE = m_pMesh->halfedgeEdge( pH );
		if( !m_pMesh->isBoundary( pE ) ) continue;
		start = i;
		break;
	}
	std::vector<CConvexCapHalfEdge*> hes;
	for( int i = 0; i < 3; i ++ )
	{
		hes.push_back( hs[(i+start)%3] );
	}

	std::vector<CConvexCapVertex*> vs;
	for( int i = 0; i < 3; i ++ )
	{
		CConvexCapHalfEdge * pH = hes[i];
		CConvexCapVertex * pV = m_pMesh->halfedgeSource( pH ); 
		vs.push_back( pV );
	}

	/* compute the length */
	std::vector<double> ls;

	for( size_t i = 0; i < hes.size(); i ++ )
	{
		CConvexCapHalfEdge * pH = hes[i];
		CConvexCapEdge     * pE = m_pMesh->halfedgeEdge( pH );		
		CConvexCapVertex * pT = m_pMesh->halfedgeTarget( pH );
		CConvexCapVertex * pS = m_pMesh->halfedgeSource( pH );
		
		double  l = pE->l();
		double hs = pS->h();
		double ht = pT->h();

		l = sqrt( l*l - (hs-ht) * (hs - ht) );
		ls.push_back( l );
	}

	std::vector<double> omega;
	for( int i = 0; i < 3; i ++ )
	{
		double d = ls[(i+1)%3]*ls[(i+1)%3] + ls[(i+2)%3] * ls[(i+2)%3] - ls[(i+0)%3] * ls[(i+0)%3];
		d = d/( 2 * ls[(i+1)%3] * ls[(i+2)%3] );
		omega.push_back( acos( d ) );
	}

	


	for( int i = 0; i < 3; i ++ )
	{
		hes[i]->omega() = omega[(i+2)%3];
	}

	CConvexCapVertex * pV = _prism.idVertex(0);
	pV->point() = CPoint(0,0,0);
	pV = _prism.idVertex(1);
	pV->point() = CPoint(ls[0],0,0);
	pV = _prism.idVertex(2);
	pV->point() = CPoint(ls[2]*cos(omega[1]),ls[2]*sin(omega[1]),0);

	pV = _prism.idVertex(2);
	CConvexCapVertex * pW = _prism.idVertex(3);

	pW->point() = pV->point() + CPoint(0,0, vs[2]->h() );

}


/*!	Compute harmonic map using direct method
*/
/*
void CConvexCap::_map()
{
	while( true )
	{
		for( CCCMesh::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
		{
			CConvexCapFace * pF = *fiter;
			bool is_boundary = false;
			for( CCCMesh::FaceVertexIterator fviter( pF ); !fviter.end(); fviter ++ )
			{
				CConvexCapVertex * pV = *fviter;
				if( m_pMesh->isBoundary( pV ) ) 
				{
					is_boundary = true;
					break;
				}
			}

			if( is_boundary )
			{
				_embed_boundary_prism( pF, m_boundary_prism );
				//m_boundary_prism.write_m("boundary_prism.m");
				__compute_prism_face_normal( m_boundary_prism );
				__compute_prism_edge_dihedral_angle( m_boundary_prism );
			}
			else
			{
				_embed_interior_prism( pF, m_interior_prism );
				//m_interior_prism.write_m("interior_prism.m");
				__compute_prism_face_normal( m_interior_prism );
				__compute_prism_edge_dihedral_angle( m_interior_prism );
			}
		}

		for( CCCMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
		{
			CConvexCapVertex * pV = *viter;

			if( m_pMesh->isBoundary( pV ) ) continue;

			pV->k() = 2 * 3.1415926535;

			for( CCCMesh::VertexInHalfedgeIterator vhiter( m_pMesh, pV); !vhiter.end(); vhiter ++ )
			{
				CConvexCapHalfEdge * pH = *vhiter;
				pV->k() -= pH->omega();
			}

			//printf("Vertex %d K %f\n", pV->id(), pV->k() );
		}

		double sum = 0;

		for( CCCMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
		{
			CConvexCapVertex * pV = *viter;

			if( m_pMesh->isBoundary( pV ) ) continue;

			pV->h() += pV->k() * 1e-4;
			sum += pV->k() * pV->k();
		}
		
		printf("%f\n", sum );

		if( sum < 1e-10 ) break;
	}
}
*/

/*!	Compute harmonic map using direct method
*/
void CConvexCap::_map()
{
	while( true )
	{

	for( CCCMesh::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		CConvexCapEdge * pE = *eiter;
		CConvexCapVertex * pV1 = m_pMesh->edgeVertex1( pE );
		CConvexCapVertex * pV2 = m_pMesh->edgeVertex2( pE );
		double h1 = pV1->h();
		double h2 = pV2->h();
		double l = pE->l();
		l = l*l - (h1-h2)*(h1-h2);
		if( l < 0 )
		{
			std::cerr << "The triangle inequality fails" << std::endl;
		}

		pE->source_l() = sqrt( l );
	}

	for( CCCMesh::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		CConvexCapFace * pF = *fiter;
		std::vector<CConvexCapHalfEdge*> hs;
		std::vector<double> ls;

		for( CCCMesh::FaceHalfedgeIterator hiter( pF ); !hiter.end(); hiter ++ )
		{
			CConvexCapHalfEdge * pH = *hiter;
			hs.push_back( pH );

			CConvexCapEdge * pE = m_pMesh->halfedgeEdge( pH );
			ls.push_back( pE->source_l() );
		}

		for( int i =  0; i < 3; i ++ )
		{
			double a = ls[(i+0)%3];
			double b = ls[(i+1)%3];
			double c = ls[(i+2)%3];

			double C = (a*a+b*b-c*c)/(2*a*b);

			if( fabs(C) > 1 ) 
			{
				std::cerr << "The triangle inequality fails" << std::endl;
			}

			C = acos(C);

			hs[i]->omega() = C;
		}


	}

	for( CCCMesh::MeshVertexIterator viter( m_pMesh); !viter.end(); viter ++ )
	{
		CConvexCapVertex * pV = *viter;
		if( m_pMesh->isBoundary( pV ) ) continue;

		pV->k() = 2 * 3.1415926535;

		for( CCCMesh::VertexInHalfedgeIterator vhiter( m_pMesh, pV ); !vhiter.end(); vhiter ++ )
		{
			CConvexCapHalfEdge * pH = *vhiter;
			pV->k() -= pH->omega();
		}
	}

	double sum = 0;
	for( CCCMesh::MeshVertexIterator viter( m_pMesh); !viter.end(); viter ++ )
	{
		CConvexCapVertex * pV = *viter;
		if( m_pMesh->isBoundary( pV ) ) continue;
		
		sum += pV->k() * pV->k();
		pV->h() += pV->k() * 1e-3;
	}

	std::cout << "Sum " << sum << std::endl;

	if( sum < 1e-9) break;
}
}
