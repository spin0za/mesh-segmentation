/**
* @file
* @author  Yuren Huang
*
* @section DESCRIPTION
*
*
* Here an isoline is represented by a vector of intersection points between 
* the real continuous isoline and the triangle edges.
*
* Essentially it's an implementation of marching triangle (a variation of the 
* marching cube algorithm), an isoline finding method. See the function 
* getIsoline() for details.
*/

#pragma once
#include "ViewerMesh/ViewerDynamicMesh.h"
#include "Operator/Operator.h"
#include <assert.h>
#include <list>
#include <time.h>
#include "cholmod.h"
#include "functions.h"
using namespace MeshLib;

typedef struct endpoints
{
	CViewerVertex* start;
	CViewerVertex* end;
}ENDS;
// note that here "L" means a vector, not list
typedef std::vector<ENDS> ENDSL;
typedef std::vector<CPoint> PL;
typedef std::vector<CViewerFace*> FL;
typedef std::vector<CViewerEdge*> EL;
typedef std::vector<CViewerVertex*> VL;
typedef std::vector<PL> MPL;
typedef std::vector<FL> MFL;
typedef std::vector<EL> MEL;
typedef std::vector<VL> MVL;
typedef std::vector<double> ARR;
typedef struct loop
{
	PL points;
	FL faces;
	EL edges;
	VL vertices;
	double length;

	void clear()
	{
		points.clear();
		faces.clear();
		edges.clear();
		vertices.clear();
		length = 0;
	}
	size_t size()
	{
		return points.size();
	}
}LOOP;
typedef struct isolines
{
	std::vector<LOOP> loops;
	double thresh;

	void push_back(LOOP l)
	{
		loops.push_back(l);
	}
	void clear()
	{
		loops.clear();
		thresh = 0;
	}
	bool empty()
	{
		return loops.empty();
	}
	size_t size()
	{
		return loops.size();
	}
	LOOP & operator[] (int i)
	{
		assert(0 <= i && i < size());
		return loops[i];
	}
}ISO;
typedef std::vector<ISO> MISO;
enum class INTERSECTION { EDGE, VERTEX, NONE };

class CIsoline
{
public:
	CIsoline();
	~CIsoline();

	void computeMeshInfo(CVDMesh *);
	void cholmodEntry(cholmod_triplet *, int, int, double, cholmod_common *);
	void generateField(CVDMesh *, ENDSL);
	ISO & getIsoline(CVDMesh *, double);
	void clearMarks(ISO, int);
	void clearMarks(ISO);
	void clearMarks();
	void clearPtrs();

private:
	ISO m_isoline;
	LOOP m_loop;
	FL m_marked_faces;
	EL m_marked_edges;
	VL m_marked_vertices;
	INTERSECTION m_crossing_type;
	int m_maxDegree;

	CIsoline(const CIsoline &);
};