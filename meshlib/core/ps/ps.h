#ifndef _PS_H_
#define _PS_H_

#include <stdio.h>
#include <stdlib.h>

#include "Riemannian/RicciFlow/RicciFlowMesh.h"
#include "Geometry/Point2.h"

namespace MeshLib
{
namespace RicciFlow
{
class CPs
{
public:
	CPs( CRFMesh * pMesh ) { m_pMesh = pMesh; m_xmin = 0; m_xmax = 1020; m_ymin = 0; m_ymax = 1020; m_eps = 1; };
	~CPs() {};
	void print( const char * filename );
	void hyperbolic_print( const char * filename );

protected:
	CRFMesh       * m_pMesh;
	double		  m_xmin; 
	double		  m_ymin; 
	double		  m_xmax; 
	double		  m_ymax; 
	int			  m_eps;

	void  _hyperbolic_euclidean( CPoint2 c, double r, CPoint2 & C, double & R ); //converting a hyperbolic circle (c,r) to a Euclidean circle (C,R)
};

} //namespace RicciFlow
} //namespace MeshLib
#endif