#ifndef RENDER_OPTIONS_WINDOW_H
#define RENDER_OPTIONS_WINDOW_H

#include <QWidget>
#include <string>
#include <iostream>
#include "MeshViewer.h"

class QGroupBox;
class QCheckBox;
class QRadioButton;

class RenderOptionsWindow : public QWidget
{
	Q_OBJECT

public:
	RenderOptionsWindow(QWidget *parent = 0);
	RenderOptionsWindow(CMeshViewer * pViewer, QWidget *parent = 0);

	void setViewer(CMeshViewer * pViewer)
	{
		m_pViewer = pViewer;
	};

private:
	QGroupBox *createFirstExclusiveGroup();
	QGroupBox *createSecondExclusiveGroup();
	QGroupBox *createNonExclusiveGroup();
	QGroupBox *createPushButtonGroup();

	QGroupBox *createTextureGroup();
	QGroupBox *createNormalGroup();
	QGroupBox *createGeometryGroup();
	QGroupBox *createPolygonGroup();
	QGroupBox *createManipulationGroup();
	QGroupBox *createShowOptions();
	QGroupBox *createPickGroup();
	QGroupBox *createBrushGroup();

private:

	QGroupBox    * texture_group;
	QRadioButton * texture_replace;
	QRadioButton * texture_modulate;
	QRadioButton * texture_none;

	QGroupBox    * normal_group;
	QRadioButton * normal_face;
	QRadioButton * normal_vertex;

	QGroupBox    * polygon_group;
	QRadioButton * polygon_fill;
	QRadioButton * polygon_line;

	QGroupBox    * geometry_group;
	QRadioButton * geometry_space;
	QRadioButton * geometry_plane;

	QGroupBox    * manipulation_group;
	QRadioButton * manipulation_obj;
	QRadioButton * manipulation_eye;

	QGroupBox    * pick_group;
	QRadioButton * pick_vertex;
	QRadioButton * pick_face;
	
	QGroupBox    * brush_group;
	QRadioButton * brush_part;
	QRadioButton * brush_patch;


	QCheckBox	 * show_edge;
	QCheckBox    * show_boundary;
	QCheckBox    * show_backface;
	QCheckBox    * show_sharp_edge;
	QCheckBox    * show_vertex_rgb;
	QCheckBox    * normal_filter;
	QCheckBox    * normal_reverse;
	QCheckBox    * show_picked;
	QCheckBox    * slow_motion;

	public slots:

	void textureGroupClicked(bool b);
	void textureclicked();
	void normalclicked();
	void polygonclicked();
	void geometryclicked();
	void manipulationclicked();
	void pickclicked();
	void brushclicked();

	void ShowEdgeClicked(bool b);
	void ShowSharpEdgeClicked(bool b);
	void ShowVertexRgbClicked(bool b);
	void ShowBoundaryClicked(bool b);
	void ShowBackFaceClicked(bool b);
	void NormalFilterClicked(bool b);
	void NormalReverseClicked(bool b);
	void ShowPickedClicked(bool b);
	void SlowMotionClicked(bool b);

private:
	CMeshViewer * m_pViewer;

};

#endif