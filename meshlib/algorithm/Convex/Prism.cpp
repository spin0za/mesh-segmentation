#include "Prism.h"

using namespace MeshLib;

void CPrism::_create_interior_prism()
{
	
	for(int i = 0; i < 6; i ++ )
	{
		CPrismVertex * pv = m_interior_prism.createVertex(i);
	}

	std::vector<CPrismVertex*> vs;

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




void CPrism::__compute_prism_face_normal( CPrMesh & _prism )
{
	for( CPrMesh::MeshFaceIterator fiter( &_prism ); !fiter.end(); fiter ++ )
	{
		CPrismFace * pF = *fiter;
		std::vector<CPrismVertex*> vs;
		for( CPrMesh::FaceVertexIterator fviter( pF ); !fviter.end(); fviter ++ )
		{
			CPrismVertex * pV = *fviter;
			vs.push_back( pV );
		}

		CPoint n = (vs[2]->point() - vs[1]->point())^(vs[0]->point() - vs[1]->point() );
		pF->normal() = n/n.norm();
	}
}

/*!
 *	Compute edge dihedral angle of prism
 */
void CPrism::__compute_prism_edge_dihedral_angle( CPrMesh & _prism )
{
	for( CPrMesh::MeshEdgeIterator eiter( &_prism ); !eiter.end(); eiter ++ )
	{
		CPrismEdge * pE = *eiter;

		CPrismHalfEdge * ph0 = _prism.edgeHalfedge(pE,0);
		CPrismHalfEdge * ph1 = _prism.edgeHalfedge(pE,1);

		CPrismFace * pF0 = _prism.halfedgeFace( ph0 );
		CPrismFace * pF1 = _prism.halfedgeFace( ph1 );
		
		CPrismVertex * pT = _prism.halfedgeTarget( ph0 );
		CPrismVertex * pS = _prism.halfedgeSource( ph0 );

		double d = pF0->normal()* pF1->normal();

		pE->angle()= 3.1415926535-acos(d);

		//std::cout << "Dihedral Angle (" << pS->id() << "," << pT->id()<<  ")" << pE->angle()/3.1415926535 << "Pi" << std::endl;
	}
}


CPrism::~CPrism()
{
};

CPrism::CPrism()
{
	_create_interior_prism();
}

void CPrism::_embed_interior_prism( std::vector<double> & ls, std::vector<double> & hs, CPrMesh & _prism )
{
	std::vector<double> bls;

	for( int i = 0;i < 3; i ++ )
	{
		double b = ls[i] * ls[i] - ( hs[(i+1)%3] - hs[(i+2)%3] ) * ( hs[(i+1)%3] - hs[(i+2)%3] );
		if( b < 0 )
		{
			std::cerr << "Error in _embed_interior_prism " << "Triangle inequality fails " << std::endl;
			return;
		}
		bls.push_back( sqrt( b ) );
	}
	
	std::vector<double> angs;

	for( int i = 0; i < 3; i ++ )
	{
		double a = bls[(i+1)%3];
		double b = bls[(i+2)%3];
		double c = bls[(i+0)%3];

		double C = ( a*a+b*b-c*c )/ ( 2*a*b );

		if( fabs(C) > 1 )
		{
			std::cerr << "Error in _embed_interior_prism " << "Cosine Law fails " << std::endl;
			return;			
		}

		C = acos( C );

		angs.push_back( C );
	}

	CPrismVertex * pV = m_interior_prism.idVertex( 0 );
	pV->point() = CPoint(0,0,0);
	
	pV = m_interior_prism.idVertex( 1 );
	pV->point() = CPoint(bls[2],0,0);

	pV = m_interior_prism.idVertex( 2 );
	pV->point() = CPoint(bls[1]*cos(angs[0]), bls[1]*sin(angs[0]),0 );

	for( int i = 0; i < 3; i ++ )
	{
		CPrismVertex * pV = m_interior_prism.idVertex( i + 0);
		CPrismVertex * pW = m_interior_prism.idVertex( i + 3);

		pW->point() = pV->point() + CPoint(0,0, hs[i]);
	}
}


void CPrism::__compute_prism_face_corner_angle( CPrMesh & _prism )
{
	
	for( CPrMesh::MeshFaceIterator fiter( &_prism ); !fiter.end(); fiter ++ )
	{
		CPrismFace * pF = *fiter;
		std::vector<CPrismHalfEdge*> hs;
		std::vector<CPrismVertex*> vs;
		
		for( CPrMesh::FaceHalfedgeIterator hiter( pF ); !hiter.end(); hiter ++ )
		{
			CPrismHalfEdge * pH = *hiter;
			hs.push_back( pH );
			vs.push_back( _prism.halfedgeTarget( pH ) );
		}
		
		CPoint e1 = vs[1]->point() - vs[0]->point();
		e1 /= e1.norm();

		CPoint n = (vs[2]->point() - vs[1]->point())^(vs[0]->point() - vs[1]->point() );
		
		CPoint e2 = n^e1;
		e2 /= e2.norm();

		std::vector<CPoint> ps;

		for( size_t i = 0; i < vs.size(); i ++ )
		{
			CPoint p = vs[i]->point();
			p -= vs[0]->point();

			CPoint q( p*e1, p*e2, 0 );
			ps.push_back( q );
		}

		for( size_t i = 0; i < ps.size(); i ++ )
		{
			CPoint q = ps[i];
			CPoint lq = ps[(i-1+ps.size())%ps.size()];
			CPoint rq = ps[(i+1+ps.size())%ps.size()];
			
			CPoint l = lq - q;
			CPoint r = rq - q;

			double ang = (l*r)/(l.norm()*r.norm());

			hs[i]->angle() = acos(ang);

			CPrismVertex * lv = vs[(i-1+vs.size())%vs.size()];
			CPrismVertex * cv = vs[i];
			CPrismVertex * rv = vs[(i+1+vs.size())%vs.size()];

			//std::cout << "Angle (" << lv->id()<< ","<< cv->id()<<","<< rv->id()<<")" << " " << hs[i]->angle()/3.1415926535 << "Pi"<< std::endl ;
		}
	}


}
