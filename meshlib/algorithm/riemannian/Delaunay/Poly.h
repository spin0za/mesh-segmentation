/*! \file Poly.h
*   \brief PLSG file - poly
*   \author David Gu
*   \date Documented on 10/14/2010
*/

#ifndef  _POLY_H_
#define  _POLY_H_

#include <map>
#include <vector>
#include <queue>

#include "Mesh/Vertex.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "Mesh/HalfEdge.h"

#include "Mesh/BaseMesh.h"
#include "Mesh/boundary.h"
#include "Mesh/iterators.h"
#include "Parser/parser.h"


namespace MeshLib
{
 /*! \brief CPolySegment class
     
	 a segment in PLSG
	 */
class CPolySegment
{
public:
	CPolySegment( int s, int t) { m_v[0] = s; m_v[1] = t;};
	int & vertex(int i) { return m_v[i]; };
protected:
	int m_v[2];
};

/*! \brief CPoly class
*
*  PLSG representation
*/
class CPoly
{
public:
	/*! CPoly default constructor
	 */
	CPoly(){};

	/*! CPoly constructor
	 * \param name input file name
	 */
	CPoly( const char * name )
	{
		std::fstream is( name, std::fstream::in );
		if( is.fail() )
		{
			fprintf(stderr,"Error in opening file %s\n", name );
			return;
		}

		char buffer[MAX_LINE];
		
		//read in # of points
		{
			is.getline(buffer, MAX_LINE );
			std::string line( buffer );
			strutil::Tokenizer stokenizer( line, " \r\n" );
			stokenizer.nextToken();
			std::string token = stokenizer.getToken();
			m_num_pts = strutil::parseString<int>(token);
		}
		for( int i = 0; i < m_num_pts; i ++ )
		{
			is.getline(buffer, MAX_LINE );
			std::string line( buffer );
			strutil::Tokenizer stokenizer( line, " \r\n" );
			 
			 //id
			 stokenizer.nextToken();
			 std::string token = stokenizer.getToken();
			 //point
			 stokenizer.nextToken();
			 token = stokenizer.getToken();
			 double x = strutil::parseString<double>(token);
			 
			 stokenizer.nextToken();
			 token = stokenizer.getToken();
			 double y = strutil::parseString<double>(token);

			 m_points.push_back(CPoint2(x,y));
		}
		{
			is.getline(buffer, MAX_LINE ); //blank line
			//read in num of segments
			is.getline(buffer, MAX_LINE ); 
			std::string line( buffer );
			strutil::Tokenizer stokenizer( line, " \r\n" );
			stokenizer.nextToken();
			std::string token = stokenizer.getToken();
			m_num_segs = strutil::parseString<int>(token);
		}
		for( int i = 0; i < m_num_segs; i ++ )
		{
			is.getline(buffer, MAX_LINE ); //blank line
			std::string line( buffer );
			strutil::Tokenizer stokenizer( line, " \r\n" );
			 //id
			 stokenizer.nextToken();
			 std::string token = stokenizer.getToken();
			 //segment id
			 stokenizer.nextToken();
			 token = stokenizer.getToken();
			 int s = strutil::parseString<int>(token);
			 
			 stokenizer.nextToken();
			 token = stokenizer.getToken();
			 int t = strutil::parseString<int>(token);

			 m_segments.push_back(CPolySegment(s,t));
		}

		is.getline(buffer, MAX_LINE ); //bland line

		//read in # of holes
		{
			is.getline(buffer, MAX_LINE );
			std::string line( buffer );
			strutil::Tokenizer stokenizer( line, " \r\n" );
			stokenizer.nextToken();

			//number of hole centers
			stokenizer.nextToken();
			std::string token = stokenizer.getToken();
			m_num_holes = strutil::parseString<int>( token );
		}
		for( int i = 0; i < m_num_holes; i ++ )
		{
			is.getline(buffer, MAX_LINE );
			std::string line( buffer );
			strutil::Tokenizer stokenizer( line, " \r\n" );
			 
			 //id
			 stokenizer.nextToken();
			 std::string token = stokenizer.getToken();
			 //point
			 stokenizer.nextToken();
			 token = stokenizer.getToken();
			 double x = strutil::parseString<double>(token);
			 
			 stokenizer.nextToken();
			 token = stokenizer.getToken();
			 double y = strutil::parseString<double>(token);

			 m_hole_centers.push_back(CPoint2(x,y));
		}

	};
	/*! CPoly destructor
	 * \param name input file name
	 */

	~CPoly(){};
	/*! vector of input points */
	std::vector<CPoint2>        & points() { return m_points; };
	/*! vector of input segment */
	std::vector<CPolySegment>   & segments() { return m_segments; };
	/*! vector of hole centers */
	std::vector<CPoint2>        & holes_centers() { return m_hole_centers; };
protected:
	/*! vector of input points */
	std::vector<CPoint2>  m_points;
	/*! vector of input segment */
	std::vector<CPolySegment> m_segments;
	/*! vector of hole centers */
	std::vector<CPoint2>  m_hole_centers;
	/*! number of points */
	int m_num_pts; 
	/*! number of segments */
	int m_num_segs;
	/*! number of holes */
	int m_num_holes;
};

}
#endif  _POLY_H_