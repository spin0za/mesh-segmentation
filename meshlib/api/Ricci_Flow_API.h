#ifndef _RICCI_FLOW_API_H_
#define _RICCI_FLOW_API_H_

#include "Operator/Operator.h"
//for topological torus
#include "Riemannian/RicciFlow/TangentialRicciFlow.h"
#include "Riemannian/RicciFlow/GenusOneTangentialRicciFlow.h"
#include "Riemannian/RicciFlow/EuclideanEmbed.h"
#include "Riemannian/RicciFlow/InversiveDistanceRicciFlow.h"
#include "Riemannian/RicciFlow/YamabeFlow.h"
#include "Riemannian/RicciFlow/VirtualInversiveDistanceRicciFlow.h"
#include "Riemannian/RicciFlow/DynamicYamabeFlow.h"
//for Extremal Lengths
#include "Riemannian/RicciFlow/TangentialRicciExtremalLength.h"
#include "Riemannian/RicciFlow/InversiveDistanceRicciExtremalLength.h"
#include "Riemannian/RicciFlow/YamabeExtremalLength.h"
#include "Riemannian/RicciFlow/VirtualInversiveDistanceRicciExtremalLength.h"
//for Riemann Mapping
#include "Riemannian/RicciFlow/TangentialRicciRiemannMapping.h"
#include "Riemannian/RicciFlow/InversiveDistanceRicciRiemannMapping.h"
#include "Riemannian/RicciFlow/YamabeRiemannMapping.h"
#include "Riemannian/RicciFlow/CurvatureFlowExponentialMap.h"
#include "Riemannian/RicciFlow/VirtualInversiveDistanceRicciPolyAnnulus.h"
//for poly annulus
#include "Riemannian/RicciFlow/TangentialRicciPolyAnnulus.h"
#include "Riemannian/RicciFlow/InversiveDistanceRicciPolyAnnulus.h"
#include "Riemannian/RicciFlow/YamabeFlowPolyAnnulus.h"
//for hyperbolic Ricci flow
#include "Riemannian/RicciFlow/TangentialHyperbolicRicciFlow.h"
#include "Riemannian/RicciFlow/InverseDistanceHyperbolicRicciFlow.h"
#include "Riemannian/RicciFlow/YamabeHyperbolicFlow.h"
#include "Riemannian/RicciFlow/VirtualInversiveDistanceRicciRiemannMapping.h"
#include "Riemannian/RicciFlow/HyperbolicEmbed.h"
/******************************************************************************************************************************
*
*	Tangential Hyperbolic RIcci Flow for Riemann Mapping
*
*******************************************************************************************************************************/
#include "Riemannian/RicciFlow/TangentialHyperbolicRicciFlowRiemannMap.h"
#include "ps/ps.h"
/******************************************************************************************************************************
*
*	QC Yamabe Flow Riemann Mapping
*
*******************************************************************************************************************************/
#include "Riemannian/RicciFlow/CurvatureFlowExponentialMap.h" 
#include "Riemannian/RicciFlow/QCInversiveDistanceRicciRiemannMapping.h"
#include "Riemannian/RicciFlow/QCYamabeRiemannMapping.h"
#include "Riemannian/RicciFlow/QCVirtualInversiveDistanceRicciRiemannMapping.h"
#include "Riemannian/RicciFlow/QCInversiveDistanceRicciRiemannMapping.h"

/*
#include "Riemannian/RicciFlow/GenusOneTangentialRicciFlow.h"
#include "Riemannian/RicciFlow/GenusOneInversiveDistanceRicciFlow.h"
#include "Riemannian/RicciFlow/InversiveDistanceRicciFlow.h"
//virtual Inversive Distance Ricci Flow
#include "Riemannian/RicciFlow/VirtualInversiveDistanceRicciFlow.h"
#include "Riemannian/RicciFlow/GraphRicciFlow.h"
#include "Topology/CoveringSpace/TorusDeckTransformation.h"
#include "Topology/CoveringSpace/Fuchs.h"
//Tangential Ricci flow for Spherical Mapping
#include "Riemannian/RicciFlow/TangentialRicciSphericalMapping.h"
//square tiling
#include "Conformal/SquareTiling/SquareTilingMesh.h"
#include "Conformal/SquareTiling/SquareTiling.h"
*/
/*!
 *	Converting from one structure to another
 */

/*******************************************************************************************************************************
 *
 *	Genus One surface
 *******************************************************************************************************************************/
//-genus_one_tangential_ricci_flow gnome.m gnome.open.m gnome.uv.m
void _genus_one_tangential_ricci_flow( const char * _closed_mesh, const char * _fundamental_domain, const char * _mesh_with_uv );
/*!
 *	Tangential Ricci Flow
 */
void _tangential_ricci_flow( const char * _closed_mesh, const char * _fundamental_domain, const char * _mesh_with_uv );
/*!
 *	Tangential Ricci Flow
 */
void _inverse_ricci_flow( const char * _closed_mesh, const char * _fundamental_domain, const char * _mesh_with_uv );
/*  Yamabe flow
 * -yamabe_flow kitten.m kitten.open.m kitten.uv.m
 */
void _Yamabe_flow( const char * _closed_mesh, const char * _fundamental_domain, const char * _mesh_with_uv );
/*
 * -virtual_inversive_ricci_flow kitten.m kitten.open.m kitten.uv.m
 */
void _virtual_inversive_ricci_flow( const char * input_mesh, const char * fundamental_domain, const char * uv_mesh );
/*  Dynamic Yamabe flow
 * -dynamic_yamabe_flow kitten.m kitten.uv.m
 */
void _dynamic_Yamabe_flow( const char * _closed_mesh, const char * _mesh_with_uv );


/*******************************************************************************************************************************
 *
 *	Extremal Length
 *
 *******************************************************************************************************************************/

/*!
 * -tangent_ricci_extremal_length sophie.remesh.m sophie.uv.m
 */
void _tangent_ricci_extremal_length( const char * _input_mesh, const char * _mesh_with_uv );

/*! 
 * -idrf_extremal_length sophie.remesh.m sophie.uv.m
 */
void _idrf_extremal_length( const char * _input_mesh, const char * _mesh_with_uv );

/*! 
 * -yamabe_extremal_length sophie.remesh.m sophie.uv.m
 */
void _Yamabe_extremal_length( const char * _input_mesh, const char * _mesh_with_uv );
/*!  
 * -virtual_extremal_length sophie.m sophie.uv.m
 */
void _virtual_inversive_distance_extremal_length( const char * input_mesh, const char * uv_mesh );


/*******************************************************************************************************************************
 *
 *	Riemann Mapper
 *
 *******************************************************************************************************************************/

/*!
 * -ricci_riemann_map sophie.m sophie.open.m sophie.uv.m
 */
void _tangent_ricci_riemann_map( const char * _closed_mesh, const char * _fundamental_domain, const char * _mesh_with_uv );
/*!
 * -idrf_riemann_map sophie.m sophie.open.m sophie.uv.m
 */
void _idrf_riemann_map(  const char * _closed_mesh, const char * _fundamental_domain, const char * _mesh_with_uv );
/*!
 * -yamabe_riemann_map sophie.m sophie.open.m sophie.uv.m
 */
void _yamabe_riemann_map( const char * _closed_mesh, const char * _fundamental_domain, const char * _mesh_with_uv );
/*!
 * -virtual_riemann_mapping sophie.m sophie.open.m sophie.uv.m
 */
void _virtual_distance_riemann_mapping( const char * input_mesh, const char * open_mesh, const char * uv_mesh );

/*******************************************************************************************************************************
 *
 *	Poly Annuli
 *
 *******************************************************************************************************************************/


/*!
 * -ricci_poly_annulus alex_hole.m alex_hole.open.m alex_hole.uv.m
 */
void _ricci_poly_annulus( const char * _input_mesh, const char * _fundamental_domain, const char * _mesh_with_uv );
/*!
 * -idrf_poly_annulus sophie.m sophie.open.m sophie.uv.m
 */
void _idrf_poly_annulus( const char * _input_mesh, const char * _fundamental_domain, const char * _mesh_with_uv );
/*! 
 * -yamabe_poly_annulus sophie.m sophie.open.m sophie.uv.m
 */
void _yamabe_poly_annulus( const char * _input_mesh, const char * _fundamental_domain, const char * _mesh_with_uv );
/*!
 * -virtual_poly_annulus alex_hole.m alex_hole.open.m alex_hole.uv.m
 */
void _virtual_distance_poly_annulus( const char * input_mesh, const char * open_mesh, const char * uv_mesh );


/*******************************************************************************************************************************
 *
 *	Hyperbolic cases
 *
 *******************************************************************************************************************************/

/*!
 * -hyper_ricci_flow eight.m eight.open.m eight.uv.m
 */
void _hyperbolic_tangential_ricci_flow( const char * closed_mesh, const char * fundamental_domain, const char * embedded_mesh );
/*!
 * -hyper_ricci_flow eight.m eight.open.m eight.uv.m
 */
void _hyperbolic_inverse_ricci_flow( const char * closed_mesh, const char * fundamental_domain, const char * embedded_mesh );
/*!
 * -hyper_yamabe_flow eight.m eight.open.m eight.uv.m
 */
void _hyperbolic_yamabe_flow( const char * closed_mesh, const char * open_mesh, const char * uv_mesh );

/******************************************************************************************************************************
*
*	Tangential Hyperbolic RIcci Flow for Riemann Mapping
*
*******************************************************************************************************************************/

/*!
 * -hyper_ricci_flow_riemann_map face.m 2.0 face.uv.m
 */
void _hyperbolic_ricci_flow_riemann_mapping( const char * input_mesh, const char * boundary_radius, const char * uv_mesh, const char * output_eps );

/******************************************************************************************************************************
*
*	Quasi-Conformal Ricci Flow Riemann Mapping
*
*******************************************************************************************************************************/
/*
 * -qc_yamabe_flow_riemann_map sophie.m sophie.open.m sophie.uv.m
 */
void _qc_yamabe_riemann_mapping( const char * input_mesh, const char * fundamental_domain, const char * uv_mesh );
/*
 * -qc_virtual_inversive_distance_riemann_map sophie.m sophie.open.m sophie.uv.m
 */
void _qc_virtual_inversive_distance_riemann_mapping( const char * input_mesh, const char * fundamental_domain, const char * uv_mesh );
/*
 * -qc_inversive_distance_riemann_map sophie.m sophie.open.m sophie.uv.m
 */
void _qc_inversive_distance_riemann_mapping( const char * input_mesh, const char * fundamental_domain, const char * uv_mesh );

#endif