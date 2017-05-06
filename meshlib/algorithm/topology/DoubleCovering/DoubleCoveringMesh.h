#ifndef  _DOUBLE_COVERING_MESH_H_
#define  _DOUBLE_COVERING_MESH_H_

#include <map>
#include <vector>
#include <queue>

#include "Mesh/Vertex.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/boundary.h"
#include "Mesh/iterators.h"
#include "Mesh/BaseMesh.h"




namespace MeshLib
{

namespace Topology
{
  
/*! \brief CDoubleCoveringMesh class
*   
*   Compute the double covering of a mesh with boundaries
*   preassumptions
*	1. no edge is a loop, namely, the two end vertices of any edge are different
*   2. for any face, at least one vertex is not on the boundary
*
*/
template<typename V, typename E, typename F, typename H>
class CDoubleCoveringMesh : public CBaseMesh<V,E,F,H>
{
public:
	typedef MeshVertexIterator<V,E,F,H>        MeshVertexIterator;
	typedef MeshFaceIterator<V,E,F,H>          MeshFaceIterator;
	typedef FaceVertexIterator<V,E,F,H>        FaceVertexIterator;
	/*!
	 *	Compute the double covering
	 */
	void DoubleCovering();
};






template<typename V, typename E, typename F, typename H>
void CDoubleCoveringMesh<V,E,F,H>::DoubleCovering()
{
	//--record the original vertices
	std::list<V *>origvertices;
	int vmax=0;
	for(MeshVertexIterator viter(this); !viter.end(); viter++)
	{
		V * v=*viter;
		origvertices.push_back(v);
		vmax = (vmax < v->id())? v->id():vmax;

		std::string line;
		std::stringstream iss(line);
		iss << "father=("<<v->id()<<")";
		
		if( v->string().length() > 0 )
		{
			v->string() += " ";
		}
		v->string() += iss.str();
	}

	//--record the original faces
	std::list<F *>origfaces;
	int fmax=0;
	for(MeshFaceIterator fiter(this); !fiter.end(); fiter++)
	{
		F * f=*fiter;
		origfaces.push_back(f);
		fmax = ( fmax < f->id() )? f->id():fmax;
	}


	std::map<int, int>corres;
	int currentid=vmax+1;
	for(std::list<V *>::iterator viter=origvertices.begin(); viter!=origvertices.end(); viter++)
	{
		V *v=*viter;
		if(!isBoundary(v))
		{
			V *newv=createVertex(currentid);
			newv->point()=v->point();
			newv->string() = v->string();
			corres[v->id()]=currentid;
			++currentid;

			std::string line;
			std::stringstream iss(line);
			iss << "father=("<<v->id()<<")";
			
			if( newv->string().length() > 0 )
			{
				newv->string() += " ";
			}
			newv->string() += iss.str();

		}
		else
			corres[v->id()]=v->id();
	}

	for(std::list<F *>::iterator fiter=origfaces.begin(); fiter!=origfaces.end(); fiter++)
	{
		F * f=*fiter;
		//---build v---//
		V * v[3];
		int i=0;
		for(FaceVertexIterator viter(f); !viter.end(); viter++)
		{
			v[i++]=idVertex(corres[(*viter)->id()]);
		}

		//reverse the vertex order//
		V *rv[3];
		rv[0]=v[0];
		rv[1]=v[2];
		rv[2]=v[1];
		F * newf = createFace(rv,++fmax);		
		newf->string() = f->string();
	}
}



typedef CDoubleCoveringMesh<CVertex, CEdge, CFace, CHalfEdge> CDCMesh;

} //namespace Topology
}
#endif