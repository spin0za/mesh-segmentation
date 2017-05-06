#include "HomotopyEmbed.h"
#include "Structure/Structure.h"
#include "Mesh/iterators.h"
#include "Eigen/Sparse"

using namespace MeshLib;

CHomotopyEmbed::CHomotopyEmbed( CHTMesh * pMesh ): m_pMesh( pMesh ), m_boundary( m_pMesh )
{
	int id = 0;

	for( CHTMesh::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		CHomotopyEdge * pE = *eiter;
		pE->idx() = id ++;
	}

	id = 0;

	for( CHTMesh::MeshVertexIterator viter( m_pMesh); !viter.end(); viter ++ )
	{
		CHomotopyVertex * pV = *viter;
		pV->idx() = id ++;
	}

};

CHomotopyEmbed::~CHomotopyEmbed()
{
};

void CHomotopyEmbed::_embed()
{

	std::vector<Eigen::Triplet<double> > A_coefficients;
	std::vector<Eigen::Triplet<double> > B_coefficients;

	
	//set the matrix A
	for( CHTMesh::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); ++ eiter )
	{
		CHomotopyEdge * pE = *eiter;
		int eid = pE->idx();
		
		CHomotopyVertex * pV1 = m_pMesh->edgeVertex1( pE );
		CHomotopyVertex * pV2 = m_pMesh->edgeVertex2( pE );

		int vid1 = pV1->idx();
		int vid2 = pV2->idx();
		
		CPoint dp = pV1->point() - pV2->point();

		for( int i = 0; i < 3; i ++ )
		{
			A_coefficients.push_back( Eigen::Triplet<double>(eid,3 * vid1 + i,  dp[i]) );
			A_coefficients.push_back( Eigen::Triplet<double>(eid,3 * vid2 + i, -dp[i]) );
		}

	}

	Eigen::SparseMatrix<double> A( m_interior_vertices, m_interior_vertices );
	A.setZero();
	A.setFromTriplets(A_coefficients.begin(), A_coefficients.end());


	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	std::cerr << "Eigen Decomposition" << std::endl;
	solver.compute(A);
	std::cerr << "Eigen Decomposition Finished" << std::endl;
	
	if( solver.info() != Eigen::Success )
	{
		std::cerr << "Waring: Eigen decomposition failed" << std::endl;
	}




};