/*! \file   DynamicYamabeFlowExtremalLength.h
 *  \brief  Dynamic Euclidean Yamabe flow for Extremal Length
 *  \author David Gu
 *  \date   documented on 10/16/2013
 *
 *	Algorithm for Dynamic Yamabe Flow for Extremal Length
 */

#ifndef _DYNAMIC_YAMABE_FLOW_EXTREMAL_LENGTH_H_
#define _DYNAMIC_YAMABE_FLOW_EXTREMAL_LENGTH_H_


#include "DynamicYamabeFlow.h"


namespace MeshLib
{

namespace RicciFlow
{

	/*! \brief CDynamicYamabeFlowExtremalLength class
	 *  
	 *  Algorithm for Dynamic Euclidean Yamabe flow for Extremal Length
	 */
	 
 template<typename M>
 class CDynamicYamabeFlowExtremalLength : public CDynamicYamabeFlow<M>
  {
  public:
	  /*! \brief CDynamicYamabeFlowExtremalLength constructor
	   *
	   *  call base class constructor 
	   */
	  CDynamicYamabeFlowExtremalLength( M * pMesh ): CDynamicYamabeFlow<M>( pMesh ){};
	
	  /*!
	   *	input marker file
	   */
	  void read_markers( const char * input_file );

  protected:


	/*!
	 *	Set the target curvature on each vertex
	 */
    virtual void    _set_target_curvature();

	/*!
	 *	markers on the boundary
	 */
	std::vector<typename M::CVertex*> m_markers;
  
 };





//set target curvature

template<typename M>
void CDynamicYamabeFlowExtremalLength<M>::_set_target_curvature()
{
	//uniformization metric
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
	  v->target_k() = 0;
  }

  for( size_t i = 0; i < m_markers.size(); i ++ )
  {
	  M::CVertex * pv = m_markers[i];
	  //this line is for rectangle
	  pv->target_k() = PI;
	  //this line is for triangle
	  //pv->target_k() = PI * 4.0/3.0;
  }
}
//set target curvature

template<typename M>
void CDynamicYamabeFlowExtremalLength<M>::read_markers( const char * input_file )
{
	std::fstream tis( input_file, std::fstream::in );

	if( tis.fail() )
	{
		std::cerr << "Error in opening file " << input_file << std::endl;
		return;
	}


	while( !tis.eof() )
	{
		int id;
		tis >> id;

		M::CVertex * pV = m_pMesh->idVertex( id );
		m_markers.push_back(pV);
	}

	
};

} //namespace RicciFlow
} //namespace MeshLib

#endif  _DYNAMIC_YAMABE_FLOW_EXTREMAL_LENGTH_H_