#include "isoline.h"
#define max(a,b) (((a) > (b)) ? (a) : (b))
/**
* Constructor that sets the initial crossing status to be no intersection.
*
* @param timemillis Number of milliseconds
*        passed since Jan 1, 1970.
*/
CIsoline::CIsoline()
{
	m_crossing_type = INTERSECTION::NONE;
}

CIsoline::~CIsoline()
{
	clearMarks();
	clearPtrs();
}

void CIsoline::computeMeshInfo(CVDMesh *m_pMesh)
{
	clock_t start, end;
	start = clock();

	for (CVDMesh::MeshFaceIterator fiter(m_pMesh); !fiter.end(); fiter++)
	{
		CViewerFace *pF = *fiter;

		CViewerHalfEdge *h0 = m_pMesh->faceMostCcwHalfEdge(pF);
		CViewerEdge *e0 = m_pMesh->halfedgeEdge(h0);
		double a = m_pMesh->edgeLength(e0);

		CViewerHalfEdge *h1 = m_pMesh->faceNextCcwHalfEdge(h0);
		CViewerEdge *e1 = m_pMesh->halfedgeEdge(h1);
		double b = m_pMesh->edgeLength(e1);

		CViewerHalfEdge *h2 = m_pMesh->faceNextCcwHalfEdge(h1);
		CViewerEdge *e2 = m_pMesh->halfedgeEdge(h2);
		double c = m_pMesh->edgeLength(e2);

		h0->angle() = acos((a * a + b * b - c * c) / (2 * a * b));
		h1->angle() = acos((b * b + c * c - a * a) / (2 * b * c));
		h2->angle() = acos((c * c + a * a - b * b) / (2 * c * a));
	}
	
	int id = 0;
	for (CVDMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++)
	{
		CViewerVertex *pV = *viter;
		pV->sid() = id++;
	}

	id = 0;
	for (CVDMesh::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); eiter++)
	{
		CViewerEdge *pE = *eiter;
		pE->sid() = id++;

		CViewerVertex *vi = m_pMesh->edgeVertex1(pE);
		CViewerVertex *vj = m_pMesh->edgeVertex2(pE);
		vi->degree()++;
		vj->degree()++;
	}

	m_maxDegree = 0;
	for (CVDMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++)
	{
		CViewerVertex *pV = *viter;
		if (pV->degree() > m_maxDegree)
		{
			m_maxDegree = pV->degree();
		}
	}

	end = clock();
	double duration = (end - start) / 100;
	std::cout << "Preprocessing time: " << duration / 10 << "s" << std::endl << std::endl;
}

void CIsoline::cholmodEntry(cholmod_triplet *T, int i, int j, double x, cholmod_common *cm)
{
	if (T->nnz >= T->nzmax && !cholmod_reallocate_triplet(2 * (T->nzmax), T, cm))
	{
		std::cerr << "Triplet reallocation failed." << std::endl;
	}

	/* append subscript */
	((int*)T->i)[T->nnz] = i;
	((int*)T->j)[T->nnz] = j;
	/* append value */
	((double*)T->x)[T->nnz] = x;

	T->nrow = max(T->nrow, (size_t)(i + 1));
	T->ncol = max(T->ncol, (size_t)(j + 1));
	T->nnz++;
}

void CIsoline::generateField(CVDMesh *m_pMesh, ENDSL stroke_ends)
{
	size_t n = m_pMesh->numVertices();
	size_t c = stroke_ends.size();
	double w = 1e3;

	// build the linear system A*phi = b
	cholmod_sparse *A, *At, *AtA;
	cholmod_triplet *A_coefficients;
	cholmod_dense *b, *Atb, *phi;
	cholmod_factor *L;

	cholmod_common common;
	cholmod_common *cm = &common;
	cholmod_start(cm);

	A_coefficients = cholmod_allocate_triplet(n + 2 * c, n, n * m_maxDegree + 2 * c, 0, CHOLMOD_REAL, cm);
	b = cholmod_zeros(n + 2 * c, 1, CHOLMOD_REAL, cm);
	Atb = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
	phi = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
	std::vector<double> d(n, 0);
	for (CVDMesh::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); eiter++)
	{
		CViewerEdge *pE = *eiter;
		CViewerHalfEdge *h0 = m_pMesh->edgeHalfedge(pE, 0);
		CViewerHalfEdge *h1 = m_pMesh->edgeHalfedge(pE, 1);
		double angle_a = m_pMesh->halfedgeNext(h0)->angle();
		double angle_b = m_pMesh->halfedgeNext(h1)->angle();
		CViewerVertex *vi = m_pMesh->edgeVertex1(pE);
		CViewerVertex *vj = m_pMesh->edgeVertex2(pE);
		// cotangent weight
		double w_ij = (1.0 / tan(angle_a) + 1.0 / tan(angle_b)) / 2.0;
		// uniform weight
		//double w_ij = 1;
		cholmodEntry(A_coefficients, vi->sid(), vj->sid(), -w_ij, cm);
		cholmodEntry(A_coefficients, vj->sid(), vi->sid(), -w_ij, cm);
		d[vi->sid()] += w_ij;
		d[vj->sid()] += w_ij;
	}
	for (int i = 0; i < n; i++)
	{
		// row sum diagonal
		cholmodEntry(A_coefficients, i, i, d[i], cm);
	}
	for (int i = 0; i < c; i++)
	{
		CViewerVertex *pi = stroke_ends[i].start;
		CViewerVertex *qi = stroke_ends[i].end;
		cholmodEntry(A_coefficients, n + i, pi->sid(), w, cm);
		cholmodEntry(A_coefficients, n + i + c, qi->sid(), w, cm);
		((double*)b->x)[n + i] = w;
	}
	A = cholmod_triplet_to_sparse(A_coefficients, A_coefficients->nnz, cm);
	// compute At*b
	At = cholmod_transpose(A, 2, cm);
	double alpha[2] = { 1, 0 };
	double beta[2] = { 0, 0 };
	cholmod_sdmult(At, 0, alpha, beta, b, Atb, cm);
	// compute At*A
	AtA = cholmod_aat(At, NULL, 0, CHOLMOD_REAL, cm);
	AtA->stype = 1;
	cholmod_sort(AtA, cm);
	
	// solve the least squares i.e. At*A*phi = At*b
	L = cholmod_analyze(AtA, cm);
	cholmod_factorize(AtA, L, cm);
	phi = cholmod_solve(CHOLMOD_A, L, Atb, cm);
	for (CVDMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++)
	{
		CViewerVertex *pV = *viter;
		pV->field() = ((double*)phi->x)[pV->sid()];
	}

	// free space
	cholmod_free_factor(&L, cm);
	cholmod_free_sparse(&A, cm);
	cholmod_free_sparse(&At, cm);
	cholmod_free_sparse(&AtA, cm);
	cholmod_free_triplet(&A_coefficients, cm);
	cholmod_free_dense(&b, cm);
	cholmod_free_dense(&Atb, cm);
	cholmod_free_dense(&phi, cm);
	
	cholmod_finish(cm);
}

ISO & CIsoline::getIsoline(CVDMesh *m_pMesh, double isoValue)
{
	clearMarks();
	clearPtrs();
	bool ending = false;
	int in_counter = 0, out_counter = 0;
	do
	{
		CViewerFace *startFace = NULL;
		CViewerEdge *startEdge = NULL, *currEdge = NULL;
		CViewerVertex *startVertex = NULL, *currVertex = NULL, *vLow = NULL, *vHigh = NULL;

		// searching for an initial iso-point
		for (CVDMesh::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); eiter++)
		{
			CViewerEdge *pE = *eiter;
			if (pE->isomark() >= 0)
			{
				continue;
			}
			m_pMesh->tellDirection(pE, vLow, vHigh);
			if (vLow->isomark() >= 0 || vHigh->isomark() >= 0)
			{
				continue;
			}
			if (vHigh->field() < isoValue || vLow->field() > isoValue)
			{
				continue;
			}
			else if (vLow->field() == isoValue)
			{
				m_crossing_type = INTERSECTION::VERTEX;
				m_loop.points.push_back(vLow->point());
				vLow->isomark() = out_counter;
				m_loop.vertices.push_back(vLow);
				startVertex = vLow;
				currVertex = vLow;
				// if low == high
				if (vHigh->field() == isoValue)
				{
					m_loop.points.push_back(vHigh->point());
					vHigh->isomark() = out_counter;
					m_loop.vertices.push_back(vHigh);
					pE->isomark() = out_counter;
					m_loop.edges.push_back(pE);
					for (CVDMesh::VertexFaceIterator fiter(vHigh); !fiter.end(); fiter++)
					{
						CViewerFace *pF = *fiter;
						pF->isomark() = out_counter;
						m_loop.faces.push_back(pF);
					}
					currVertex = vHigh;
					startEdge = pE;
				}
				break;
			}
			else if (vHigh->field() == isoValue)
			{
				m_crossing_type = INTERSECTION::VERTEX;
				m_loop.points.push_back(vHigh->point());
				vHigh->isomark() = out_counter;
				m_loop.vertices.push_back(vHigh);
				startVertex = vHigh;
				currVertex = vHigh;
				break;
			}
			else
			{
				m_crossing_type = INTERSECTION::EDGE;
				double ratio = (isoValue - vLow->field()) / (vHigh->field() - vLow->field());
				CPoint isoPoint = vLow->point() * (1 - ratio) + vHigh->point() * ratio;
				pE->iso_cross() = isoPoint;
				m_loop.points.push_back(isoPoint);
				pE->isomark() = out_counter;
				m_loop.edges.push_back(pE);
				startEdge = pE;
				currEdge = pE;
				break;
			}
		}

		if (m_crossing_type == INTERSECTION::NONE)
		{
			//std::cerr << "No more isoline exists for the value " << isoValue << std::endl;
			break;
		}

		// completing the loop
		do
		{
			switch (m_crossing_type)
			{
			case INTERSECTION::VERTEX:
				m_crossing_type = INTERSECTION::NONE;
				// isoline passing through two vertices
				for (CVDMesh::VertexVertexIterator witer(currVertex); !witer.end(); witer++)
				{
					CViewerVertex *pW = *witer;
					CViewerEdge *pE = m_pMesh->vertexEdge(pW, currVertex);
					// reaching the end
					if (pW == startVertex && pE != startEdge)
					{
						for (CVDMesh::VertexFaceIterator fiter(pW); !fiter.end(); fiter++)
						{
							CViewerFace *pF = *fiter;
							pF->isomark() = out_counter;
							m_loop.faces.push_back(pF);
						}
						ending = true;
						m_loop.points.push_back(pW->point());
						break;
					}
					if (pW->isomark() >= 0)
					{
						continue;
					}
					if (pW->field() == isoValue)
					{
						m_crossing_type = INTERSECTION::VERTEX;
						m_loop.points.push_back(pW->point());
						pW->isomark() = out_counter;
						m_loop.vertices.push_back(pW);
						pE->isomark() = out_counter;
						m_loop.edges.push_back(pE);
						for (CVDMesh::VertexFaceIterator fiter(pW); !fiter.end(); fiter++)
						{
							CViewerFace *pF = *fiter;
							pF->isomark() = out_counter;
							m_loop.faces.push_back(pF);
						}
						if (currVertex == startVertex)
						{
							startEdge = pE;
						}
						currVertex = pW;
						break;
					}
				}
				if (m_crossing_type != INTERSECTION::NONE || ending)
				{
					break;
				}

				// isoline passing through a vertex and an edge
				for (CVDMesh::VertexOutHalfedgeIterator hiter(m_pMesh, currVertex); !hiter.end(); hiter++)
				{
					CViewerHalfEdge *pH = m_pMesh->halfedgeNext(*hiter);
					CViewerEdge *pE = m_pMesh->halfedgeEdge(pH);
					CViewerFace *pF = m_pMesh->vEdgeFace(currVertex, pE);
					// reaching the end
					if (pE == startEdge && pF != startFace)
					{
						pF->isomark() = out_counter;
						m_loop.faces.push_back(pF);
						ending = true;
						m_loop.points.push_back(m_loop.points.front());
						break;
					}
					if (pE->isomark() >= 0)
					{
						continue;
					}
					m_pMesh->tellDirection(pE, vLow, vHigh);
					if (vHigh->field() < isoValue || vLow->field() > isoValue)
					{
						continue;
					}
					else
					{
						m_crossing_type = INTERSECTION::EDGE;
						double ratio = (isoValue - vLow->field()) / (vHigh->field() - vLow->field());
						CPoint isoPoint = vLow->point() * (1 - ratio) + vHigh->point() * ratio;
						pE->iso_cross() = isoPoint;
						m_loop.points.push_back(isoPoint);
						pE->isomark() = out_counter;
						m_loop.edges.push_back(pE);
						pF->isomark() = out_counter;
						m_loop.faces.push_back(pF);
						if (currVertex == startVertex)
						{
							startFace = pF;
						}
						currEdge = pE;
						break;
					}
				}
				break;

			case INTERSECTION::EDGE:
				m_crossing_type = INTERSECTION::NONE;
				// isoline passing through a vertex and an edge
				CViewerVertex *pW = m_pMesh->eNextVertex(currEdge, m_pMesh->edgeFace1(currEdge));
				if (!m_loop.faces.empty())
				{
					pW = m_pMesh->eNextVertex(currEdge, m_loop.faces.back());
				}
				CViewerFace *pF = m_pMesh->vEdgeFace(pW, currEdge);
				// reaching the end
				if (pW == startVertex && pF != startFace)
				{
					for (CVDMesh::VertexFaceIterator fiter(pW); !fiter.end(); fiter++)
					{
						CViewerFace *pG = *fiter;
						pG->isomark() = out_counter;
						m_loop.faces.push_back(pG);
					}
					ending = true;
					m_loop.points.push_back(pW->point());
					break;
				}
				if (pW->field() == isoValue)
				{
					m_crossing_type = INTERSECTION::VERTEX;
					m_loop.points.push_back(pW->point());
					pW->isomark() = out_counter;
					m_loop.vertices.push_back(pW);
					for (CVDMesh::VertexFaceIterator fiter(pW); !fiter.end(); fiter++)
					{
						CViewerFace *pG = *fiter;
						pG->isomark() = out_counter;
						m_loop.faces.push_back(pG);
					}
					if (currEdge == startEdge)
					{
						startFace = pF;
					}
					currVertex = pW;
					break;
				}
				// isoline passing through two edges
				else if (pW->field() < isoValue)
				{
					m_crossing_type = INTERSECTION::EDGE;
					CViewerEdge *pE = m_pMesh->vertexEdge(pW, vHigh);
					// reaching the end
					if (pE == startEdge && pF != startFace)
					{
						pF->isomark() = out_counter;
						m_loop.faces.push_back(pF);
						ending = true;
						m_loop.points.push_back(m_loop.points.front());
						break;
					}
					double ratio = (isoValue - pW->field()) / (vHigh->field() - pW->field());
					CPoint isoPoint = pW->point() * (1 - ratio) + vHigh->point() * ratio;
					pE->iso_cross() = isoPoint;
					m_loop.points.push_back(isoPoint);
					pE->isomark() = out_counter;
					m_loop.edges.push_back(pE);
					pF->isomark() = out_counter;
					m_loop.faces.push_back(pF);
					if (currEdge == startEdge)
					{
						startFace = pF;
					}
					currEdge = pE;
					m_pMesh->tellDirection(currEdge, vLow, vHigh);
				}
				else
				{
					m_crossing_type = INTERSECTION::EDGE;
					CViewerEdge *pE = m_pMesh->vertexEdge(pW, vLow);
					// reaching the end
					if (pE == startEdge && pF != startFace)
					{
						pF->isomark() = out_counter;
						m_loop.faces.push_back(pF);
						ending = true;
						m_loop.points.push_back(m_loop.points.front());
						break;
					}
					double ratio = (isoValue - vLow->field()) / (pW->field() - vLow->field());
					CPoint isoPoint = vLow->point() * (1 - ratio) + pW->point() * ratio;
					pE->iso_cross() = isoPoint;
					m_loop.points.push_back(isoPoint);
					pE->isomark() = out_counter;
					m_loop.edges.push_back(pE);
					pF->isomark() = out_counter;
					m_loop.faces.push_back(pF);
					if (currEdge == startEdge)
					{
						startFace = pF;
					}
					currEdge = pE;
					m_pMesh->tellDirection(currEdge, vLow, vHigh);
				}
				break;
			}
			in_counter++;
		} while (!ending);

		for (int i = 0; i < m_loop.size(); i++)
		{
			CPoint p0 = m_loop.points[i];
			CPoint p1 = m_loop.points[(i + 1) % m_loop.size()];
			m_loop.length += distance(p0, p1);
		}
		m_isoline.push_back(m_loop);
		m_marked_faces.insert(m_marked_faces.end(), m_loop.faces.begin(), m_loop.faces.end());
		m_marked_edges.insert(m_marked_edges.end(), m_loop.edges.begin(), m_loop.edges.end());
		m_marked_vertices.insert(m_marked_vertices.end(), m_loop.vertices.begin(), m_loop.vertices.end());
		m_loop.clear();
		m_crossing_type = INTERSECTION::NONE;
		ending = false;
		in_counter = 0;
		out_counter++;

	} while (1);

	return m_isoline;
}

void CIsoline::clearMarks(ISO isoline, int i)
{
	for (int j = 0; j < isoline[i].faces.size(); j++)
	{
		isoline[i].faces[j]->isomark() = -1;
	}
	for (int j = 0; j < isoline[i].edges.size(); j++)
	{
		isoline[i].edges[j]->isomark() = -1;
		isoline[i].edges[j]->iso_cross() = CPoint(0, 0, 0);
	}
	for (int j = 0; j < isoline[i].vertices.size(); j++)
	{
		isoline[i].vertices[j]->isomark() = -1;
	}
}

void CIsoline::clearMarks(ISO isoline)
{
	for (int i = 0; i < isoline.size(); i++)
	{
		clearMarks(isoline, i);
	}
}

void CIsoline::clearMarks()
{
	for (int i = 0; i < m_marked_faces.size(); i++)
	{
		m_marked_faces[i]->isomark() = -1;
	}
	m_marked_faces.clear();
	for (int i = 0; i < m_marked_edges.size(); i++)
	{
		m_marked_edges[i]->isomark() = -1;
		m_marked_edges[i]->iso_cross() = CPoint(0, 0, 0);
	}
	m_marked_edges.clear();
	for (int i = 0; i < m_marked_vertices.size(); i++)
	{
		m_marked_vertices[i]->isomark() = -1;
	}
	m_marked_vertices.clear();
}

void CIsoline::clearPtrs()
{
	m_loop.clear();
	m_isoline.clear();
}