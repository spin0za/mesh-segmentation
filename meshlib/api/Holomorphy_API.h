#ifndef _HOLOMORPHY_API_H_

/*!
 *	Trait IO 
 */
#include "parser/traits_io.h"

/*************************************************************************************************************************
*
*	Conformal Geometric methods
*
*************************************************************************************************************************/


/*!	Harmonic Mapping
 */
#include "Conformal/HarmonicMapper/DiskHarmonicMapper.h" //Harmonic Mapping
/*!	Spherical Harmoinc Mapping
 */
#include "Conformal/SphericalHarmonicMap/SphericalHarmonicMap.h" //Harmonic Mapping
/*! Extremal length , harmonic form
*/
#include "Conformal/ExtremalLength/ExtremalLength.h"
/*! Combinatorial Extremal length , harmonic form
*/
#include "Conformal/ExtremalLength/CombinatorialExtremalLength.h"
/*! SlitMap
*/
#include "Conformal/SlitMap/SlitMap.h"
/*!	Harmonic Exact Form
 */
#include "Conformal/ExactForm/HarmonicExactForm.h" //harmonic exact form
/*! Holomorphic 1-form
 */
#include "Conformal/HolomorphicForm/HolomorphicForm.h" //holomorphic 1-form
/*!	Compute the Exponential Map : \f$ z\to e^{z}\f$
 */
#include "Conformal/PolarMap/PolarMap.h" //Compute the exponential map z \to e^{z}

/*!	Heat Flow
 */
#include "Conformal/ClosedForm/HeatFlow.h" 
/*!	Combinatorial Heat Flow
 */
#include "Conformal/ClosedForm/CombinatorialHeatFlow.h" 

/*!	QC Harmonic Exact Form
 */
#include "Conformal/ExactForm/QCHarmonicExactForm.h" //harmonic exact form
/*! QC Holomorphic 1-form
 */
#include "Conformal/HolomorphicForm/QCHolomorphicForm.h" //holomorphic 1-form
/*!	Heat Flow under quasi-conformal structure
 */
#include "Conformal/ClosedForm/QCHeatFlow.h" 


/************************************************************************************************************************************
*
*	Conformal Geometric Algorithms based on Diagonal Ratio
*
************************************************************************************************************************************/

/*!	Holomorphic Form based on diagonal ratio
 */
#include "Conformal/DiagonalRatio/DiagonalRatioHolomorphicForm.h"
/*!	Conformal mapping based on diagonal ratio
 */
#include "Conformal/DiagonalRatio/DiagonalRatioMapping.h"


/*!	Holomorphic quadratic form for polygon
 */
#include "Conformal/HolomorphicQuadraticForm/HolomorphicQuadraticForm.h"

/************************************************************************************************************************************
*
*	Conformal Geometric Algorithms
*
************************************************************************************************************************************/


/*! 
 * Compute harmonic mapping from topological disk to the planar unit disk
 */
void _harmonic_map( const char * _input, const char * _output );

/*! 
 * Compute spherical harmonic mapping from a topological sphere to the Euclidean unit sphere
 */
void _spherical_harmonic_map( const char * _input, const char * _output );

/*!	compute the harmonic form for extremal length, need to compare to holomorphic form, and integrate
 *
 */
void _extremal_length( const char * _input, const char * _output );

/*!	compute the harmonic form for extremal length using combinatorial Laplacian
 *
 */
void _combinatorial_extremal_length( const char * _input, const char * _output );

/*!	compute the slit map for topological multi-hole anuli to a planar annulus with concentric circular slits
 *
 */
void _slit_map( int argc, char * argv[] );

/*! compute the basis for all the exact harmonic forms on a surface
 *
 */
void _exact_form( const char * _input, const char * _output );
/*!	compute the basis of all holomorphic 1-forms on a surface
 *
 */
void _holomorphic_form( int argc,  char * argv[] );
/*!	Integrate a holomorphic 1-form on a simply connected surface
 *
 */
void _integration( const char * _form_input, const char * _domain_input, const char *_output );
/*! Compute the exponential map
 *
 */
void _polar_map( const char * _closed_mesh, const char * _open_mesh, const char * _output );
/*!	diffuse the closed one form basis to a harmonic one form
 *
 */
void _diffuse( const char * _closed_1_form, const char * _fundamental_domain, const char * _harmonic_1_form, const char * _harmonic_form_integrated );
/*!	diffuse the closed one form basis to a harmonic one form for square tiling purpose
 *
 */
void _combinatorial_diffuse( const char * _closed_1_form, const char * _fundamental_domain, const char * _harmonic_1_form, const char * _harmonic_form_integrated );

/*! compute the basis for all the exact harmonic forms on a surface
 *
 */
void _qc_exact_form( const char * _input, const char * _output );
/*!	diffuse the closed one form basis to a quasi-conformal harmonic one form
 *
 */

void _qc_diffuse( const char * _closed_1_form, const char * _fundamental_domain, const char * _harmonic_1_form, const char * _ingetrated );
/*!	compute the basis of all holomorphic 1-forms on a surface
 *
 */
void _qc_holomorphic_form( int argc, char * argv[] );
/*!	Linear combination of holomorphic 1-forms
*
*/

void _linear_combination( int argc, char * argv[] );



/**********************************************************************************************************************************************
*
*	Conformal Mapping Based on Diagnal Ratio
*	
*
**********************************************************************************************************************************************/

/*!	compute the basis of all holomorphic 1-forms based on diagonal ratio
 *
 */

void _diagonal_ratio_holomorphic_form( int argc, char * argv[] );
/*!	compute the conformal mapping based on diagonal ratio
 *
 */

void _diagonal_ratio_free_boundary_mapping( const char * _input_mesh, const char * _output_mesh );

/*!
 *	Holomorphic quadratic differential forms for polygonal shapes
 */
void _holomorphic_quadratic_form( int argc,  char * argv[] );

#endif _HOLOMORPHY_API_H_
