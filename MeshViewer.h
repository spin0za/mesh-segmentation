#ifndef _MESH_VIEWER_H
#define _MESH_VIEWER_H

#include <QGLWidget>
#include <QApplication>
#include <QDesktopWidget>
#include <QMouseEvent>
#include <QFileDialog>
#include <QMessageBox>
#include <QVector3D>
#include <QColorDialog>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>
#include <atlstr.h>
#include "functions.h"
#include "OpenGLHeader.h"
#include "ViewerMesh/ViewerDynamicMesh.h"
#include "viewer/Arcball.h"                           /*  Arc Ball  Interface         */
#include "Operator/Operator.h"
#include "parser/traits_io.h"
#include "bmp/RgbImage.h"
#include "isoline.h"

using namespace MeshLib;

enum class BRUSH_MODE { PART, PATCH };
enum class DRAW_MODE { NONE, POINTS, WIREFRAME, FLATLINES, FLAT, SMOOTH, TEXTURE };
enum class TEXTURE_MODE { NONE, REPLACE, MODULATE };
enum class NORMAL_MODE { VERTEX, FACE };
enum class POLYGON_MODE { FILL, LINE };
enum class ORIENTATION_MODE { CCW, CLW };
enum class GEOMETRY_MODE { GEOMETRY, PARAMETER };
enum class MANIPULATION_MODE { EYE, OBJECT };
enum class PICK_MODE { VERTEX, FACE };

#define EPSILON 1e-4
#define EPSILON_COS_THETA_FOR_SELECT_VERTEX 1-(1e-4)

typedef struct endfaces
{
	CViewerFace* start;
	CViewerFace* end;
}ENDSF;
typedef std::vector<ENDSF> ENDSFL;
typedef std::vector<QPoint> QPL;
typedef std::vector<QPoint> QPL;
typedef std::vector<QPL> MQPL;
typedef struct strokes
{
	QPL curr;
	MQPL all;
	// for part use
	ENDS curr_ends;
	ENDSL all_ends;
	FL curr_faces;
	MFL all_faces;
	// for patch use
	ENDSF curr_end_faces;
	ENDSFL all_end_faces;

	void part_push_back()
	{
		all.push_back(curr);
		all_ends.push_back(curr_ends);
		all_faces.push_back(curr_faces);
	}
	void part_pop_back()
	{
		all.pop_back();
		all_ends.pop_back();
		all_faces.pop_back();
	}
	void part_clear()
	{
		all.clear();
		all_ends.clear();
		all_faces.clear();
	}
	void patch_push_back()
	{
		all.push_back(curr);
		all_end_faces.push_back(curr_end_faces);
	}
	void patch_pop_back()
	{
		all.pop_back();
		all_end_faces.pop_back();
	}
	void clearcurr()
	{
		curr.clear();
		curr_faces.clear();
	}
	void clearall()
	{
		all.clear();
		all_ends.clear();
		all_faces.clear();
		all_end_faces.clear();
	}
	void clear()
	{
		clearcurr();
		clearall();
	}
	bool empty()
	{
		return all.empty();
	}
}STR;


class CMeshViewer : public QGLWidget
{
	Q_OBJECT
public:
	CMeshViewer(QWidget *_parent = 0);
	~CMeshViewer();

	void loadMeshFile(const char *, std::string sExt);
	void loadTextureFile(const char *, std::string sExt);
	void set_texture_mode(TEXTURE_MODE mode);

	RgbImage * get_texture_image()
	{
		return m_pTexture;
	}

	TEXTURE_MODE get_texture_mode()
	{
		return m_texture_mode;
	}

	void set_polygon_mode(POLYGON_MODE mode);
	POLYGON_MODE get_polygon_mode()
	{
		return m_polygon_mode;
	}
	void set_orientation_mode(ORIENTATION_MODE mode);

	void set_geometry_mode(GEOMETRY_MODE mode)
	{
		m_geometry_mode = mode;
		updateGL();
		makeCurrent();
	};
	GEOMETRY_MODE get_geometry_mode()
	{
		return m_geometry_mode;
	}

	void set_normal_mode(NORMAL_MODE mode)
	{
		m_normal_mode = mode;
		updateGL();
		makeCurrent();
	};
	NORMAL_MODE get_normal_mode()
	{
		return m_normal_mode;
	}

	void set_manipulation_mode(MANIPULATION_MODE mode)
	{
		m_manipulation_mode = mode;
		updateGL();
		makeCurrent();
	};
	MANIPULATION_MODE get_manipulation_mode()
	{
		return m_manipulation_mode;
	}
	void set_edge_flag(bool flag);
	bool get_edge_flag()
	{
		return m_show_edge;
	}

	void set_boundary_flag(bool flag);
	bool get_boundary_flag()
	{
		return m_show_boundary;
	};
	void set_backface_flag(bool flag);
	bool get_backface_flag()
	{
		return m_show_backface;
	};

	void set_filter_normal_flag(bool flag);
	bool get_filter_normal_flag()
	{
		return m_doNormalFilter;
	};

	void set_reverse_normal_flag(bool flag);
	bool get_reverse_normal_flag()
	{
		return m_reverse_normal;
	};

	void set_sharp_edge_flag(bool flag);
	bool get_sharp_edge_flag()
	{
		return m_show_sharp_edge;
	};

	void set_vertex_rgb_flag(bool flag);
	bool get_vertex_rgb_flag()
	{
		return m_show_vertex_rgb;
	};
	void set_face_rgb_flag(bool flag);
	bool get_face_rgb_flag()
	{
		return m_show_face_rgb;
	};
	void set_pick_mode(PICK_MODE mode);
	PICK_MODE get_pick_mode()
	{
		return m_pick_mode;
	};
	void set_brush_mode(BRUSH_MODE mode);
	BRUSH_MODE get_brush_mode()
	{
		return m_brush_mode;
	};
	void set_show_cut_flag(bool flag);

	void set_picked_flag(bool flag);
	bool get_picked_flag()
	{
		return m_show_picked;
	};

	void set_slow_motion_flag(bool flag);
	bool get_slow_motion_flag()
	{
		return m_slowMotion;
	};
	void set_update_mode(bool flag);
	bool get_update_mode()
	{
		return m_updateMode;
	};

	void loadFromMainWindow(std::string fileName, std::string fExt);

	public slots:

	void openMesh();
	void openTexture();
	void showPoints();
	void showWireframe();
	void showSmooth();
	void showFlat();
	void showFlatlines();
	void showBoundary();
	void showTexture();
	void saveImage();
	void enterSelectionMode();
	void quitSelectionMode();
	void enterBrushMode();
	void quitBrushMode();
	void showField();
	void showIsolines();
	void showCutResult();

private:

	void resizeGL(int, int);
	void paintGL();
	void initializeGL();
	void initializeTexture();

	void updateProjectionMatrix();

	float fovy()  const { return 40.00f; };
	float zNear() const { return 0.001f; };
	float zFar()  const { return 100.0f; };

private:

	void drawWireframe();
	void drawMesh();
	void drawVertices();
	void drawEdges();
	void drawBoundary(CVDMesh * pMesh);
	void drawBackface(CVDMesh * pMesh);
	void drawSharpEdge(CVDMesh * pMesh);
	void drawTexture();

	void setupObject();
	void setupEye();
	void setupLight();

	void snapShot(std::string fileName, std::string fileExt);

private:
	/* rotation quaternion and translation vector for the object */
	CQrot       m_ObjRot;
	CPoint      m_ObjTrans;
	CQrot       m_EyeRot;
	CPoint      m_EyeTrans;

	CVDMesh  *  m_pMesh;

	/* mesh name */

	QString     m_filename;
	std::string sFilename;

	RgbImage *  m_pTexture;

private:

	void mousePressEvent(QMouseEvent * mouseEvent);
	void mouseMoveEvent(QMouseEvent * mouseEvent);
	void mouseReleaseEvent(QMouseEvent * mouseEvent);
	void wheelEvent(QWheelEvent * mouseEvent);

	int  gButton;
	QPoint latestMousePos;
	CArcball arcball;


	bool      m_show_boundary;
	bool      m_show_edge;
	bool      m_show_vertex_rgb;
	bool      m_show_face_rgb;
	bool      m_show_backface;
	bool      m_show_sharp_edge;
	bool      m_reverse_normal;
	bool      m_show_picked;

	GLuint m_texName;

	BRUSH_MODE	  m_brush_mode;
	DRAW_MODE     m_display_mode;
	NORMAL_MODE   m_normal_mode;
	TEXTURE_MODE  m_texture_mode;
	GEOMETRY_MODE m_geometry_mode;
	POLYGON_MODE  m_polygon_mode;
	MANIPULATION_MODE m_manipulation_mode;
	PICK_MODE         m_pick_mode;

private:
	bool m_updateMode;
	bool m_slowMotion;
	bool m_showCut;
	bool m_isBrushMode;
	bool m_isSelectionMode;
	bool m_doNormalFilter;
	//for picking
	void pickRay(double x, double y);
	CViewerVertex * pickVertex();
	CViewerFace *   pickFace();
	void drawSelectedVertcies();

	// cutting related
	// part
	void computeField();
	int computeIsolines();
	int cutMesh();
	void filterIsolines();
	void dyeNeighborhood(CViewerVertex*);
	FL findPath(CViewerFace*, CViewerFace*);
	// patch
	void regionGrow();
	double cost(CViewerFace*, CViewerFace*);
	double costNF(CViewerFace*, CViewerFace*);
	void normalFilter();
	double GaussianWeight(CViewerFace*, CViewerFace*);

	QVector3D m_ray[2];

	FL m_picked_faces;
	VL m_picked_vertices;

	CIsoline *m_isoline_finder;
	ISO m_isoline;
	MISO m_isoline_sequence;
	ARR m_isoline_metrics;
	MFL m_boundary_faces;
	STR m_strokes;
	int m_selection;
	int m_part_count;
	int m_patch_count;
	int m_isoline_set_size;
	int m_cut_mesh_status;

	PL m_palette;
	CPoint m_default_rgb;
};


#endif