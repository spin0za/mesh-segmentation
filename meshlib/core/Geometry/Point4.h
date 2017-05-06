/*! \file Point4.h
*   \brief Point4
*   \author David Gu
*   \date   documented on 02/17/2011
*
*   four dimensional Point class
*/
#ifndef  _POINT_4_H_
#define  _POINT_4_H_

#include <map>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

namespace MeshLib
{

	/******************************************************************************************************
			
		CPoint4 in R4 class

	*******************************************************************************************************/

	class CPoint4
	{
	public:
		CPoint4(){};
		CPoint4( const CPoint4 & pt ) { for( int i = 0; i < 4; i ++ ) m_v[i] = pt.m_v[i]; };
		CPoint4( double x, double y, double z, double w) { m_v[0] = x; m_v[1] = y; m_v[2] = z; m_v[3] = w; };
		~CPoint4(){};
		
		double & operator[]( int i ) { return m_v[i]; };
		
		CPoint4 operator-( CPoint4 & p )
		{
			return CPoint4( m_v[0]-p[0], m_v[1]-p[1], m_v[2] - p[2], m_v[3]-p[3] );
		};

		CPoint4 operator+( CPoint4 & p )
		{
			return CPoint4( m_v[0]+p[0], m_v[1]+p[1], m_v[2] + p[2], m_v[3]+p[3] );
		};

		CPoint4 operator+( double s )
		{
			return CPoint4( m_v[0] * s, m_v[1] * s, m_v[2] * s, m_v[3] * s );
		};

		CPoint4 operator/( double s )
		{
			return CPoint4( m_v[0] / s, m_v[1] / s, m_v[2] / s, m_v[3] / s );
		};

	protected:
		double m_v[4];
	};
	

}

#endif