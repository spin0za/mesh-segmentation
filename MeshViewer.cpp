#include "MeshViewer.h"
#define INFI 1e8
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define CYAN	CPoint(  0, 127, 255) / 255
#define GOLDEN	CPoint(255, 255,   0) / 255
#define GREEN	CPoint(  0, 255,   0) / 255
#define RED		CPoint(255,   0,   0) / 255
#define	INDIGO  CPoint( 75,   0, 130) / 255
#define ORANGE  CPoint(255, 165,   0) / 255
#define VIOLET  CPoint(127,   0, 255) / 255

CMeshViewer::CMeshViewer(QWidget *_parent) : QGLWidget(_parent)
{
	setAttribute(Qt::WA_NoSystemBackground, true);
	setFocusPolicy(Qt::StrongFocus);

	m_isoline_set_size = 15;
	m_selection = -1;
	m_part_count = 1;
	m_patch_count = 1;
	m_showCut = false;
	m_slowMotion = false;
	m_updateMode = false;
	m_brush_mode = BRUSH_MODE::PART;
	m_isoline_finder = nullptr;
	// m_default_rgb = CPoint(229, 162, 141) / 255;
	m_default_rgb = CPoint(204, 204, 204) / 255;
	m_isBrushMode = false;
	m_isSelectionMode = false;
	m_show_face_rgb = false;
	m_doNormalFilter = true;

	m_palette.push_back(CYAN);
	m_palette.push_back(GOLDEN);
	m_palette.push_back(GREEN);
	m_palette.push_back(RED);
	m_palette.push_back(INDIGO);
	m_palette.push_back(ORANGE);
	m_palette.push_back(VIOLET);

	m_pMesh = nullptr;
	m_pTexture = nullptr;
	m_display_mode = DRAW_MODE::SMOOTH;
	m_texture_mode = TEXTURE_MODE::NONE;
	m_geometry_mode = GEOMETRY_MODE::GEOMETRY;
	m_normal_mode = NORMAL_MODE::VERTEX;
	m_polygon_mode = POLYGON_MODE::FILL;
	m_manipulation_mode = MANIPULATION_MODE::OBJECT;
	m_pick_mode = PICK_MODE::VERTEX;

	m_show_boundary = false;
	m_show_edge = false;
	m_show_sharp_edge = false;
	m_show_backface = true;
	m_reverse_normal = false;
	m_show_vertex_rgb = true;
	m_show_picked = true;

	m_EyeTrans = CPoint(0, 0, 5);
}


CMeshViewer::~CMeshViewer()
{
	if (m_isoline_finder)
	{
		delete m_isoline_finder;
		m_isoline_finder = nullptr;
	}
	if (m_pMesh)
	{
		delete m_pMesh;
		m_pMesh = nullptr;
	}
}

void CMeshViewer::drawMesh()
{
	if (!m_pMesh) return;

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0, 1.0);

	switch (m_texture_mode)
	{
	case TEXTURE_MODE::REPLACE:
	case TEXTURE_MODE::MODULATE:
		glBindTexture(GL_TEXTURE_2D, m_texName);
		break;
	case TEXTURE_MODE::NONE:
		break;
	}
	
	// display the cutting boundary for part brush
	if (!m_boundary_faces.empty() && m_showCut)
	{
		CViewerEdge *bdCrosses[2], *e0, *e1, *e2;
		CViewerVertex *v[3], *v0, *v1, *v2;
		CPoint n;
		bool oriented;
		for (int i = 0; i < m_boundary_faces.size(); i++)
		{
			for (int j = 0; j < m_boundary_faces[i].size(); j++)
			{
				CViewerFace *pF = m_boundary_faces[i][j];
				int edgeIsoCount = 0;
				for (CVDMesh::FaceHalfedgeIterator hiter(pF); !hiter.end(); hiter++)
				{
					CViewerEdge *pE = m_pMesh->halfedgeEdge(*hiter);
					if (pE->bd_cross() != CPoint(0, 0, 0))
					{
						bdCrosses[edgeIsoCount] = pE;
						edgeIsoCount++;
					}
				}

				switch (m_display_mode)
				{
				case DRAW_MODE::FLAT:
				case DRAW_MODE::FLATLINES:
					m_normal_mode = NORMAL_MODE::FACE;
					break;
				case DRAW_MODE::SMOOTH:
					m_normal_mode = NORMAL_MODE::VERTEX;
					break;
				}

				int a = m_reverse_normal ? -1 : 1;
				double ratio = 0;
				switch (edgeIsoCount)
				{
				case 2:
					e0 = bdCrosses[0];
					e1 = bdCrosses[1];
					oriented = m_pMesh->orientEdges(e0, e1);
					assert(oriented);
					v2 = m_pMesh->eCommonVertex(e0, e1);
					v1 = m_pMesh->eOtherVertex(e0, v2);
					v0 = m_pMesh->eOtherVertex(e1, v2);

					glColor3f(v2->rgb()[0], v2->rgb()[1], v2->rgb()[2]);
					glBegin(GL_TRIANGLES);
					switch (m_normal_mode)
					{
					case NORMAL_MODE::VERTEX:
						n = v2->normal();
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(v2->point()[0], v2->point()[1], v2->point()[2]);

						ratio = (e1->bd_cross() - v2->point()).norm() / (v0->point() - v2->point()).norm();
						n = v2->normal() * ratio + v0->normal() * (1 - ratio);
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(e1->bd_cross()[0], e1->bd_cross()[1], e1->bd_cross()[2]);

						ratio = (e0->bd_cross() - v2->point()).norm() / (v1->point() - v2->point()).norm();
						n = v2->normal() * ratio + v1->normal() * (1 - ratio);
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(e0->bd_cross()[0], e0->bd_cross()[1], e0->bd_cross()[2]);
						break;
					case NORMAL_MODE::FACE:
						n = pF->normal();
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(v2->point()[0], v2->point()[1], v2->point()[2]);
						glVertex3d(e1->bd_cross()[0], e1->bd_cross()[1], e1->bd_cross()[2]);
						glVertex3d(e0->bd_cross()[0], e0->bd_cross()[1], e0->bd_cross()[2]);
						break;
					};
					glEnd();

					glColor3f(v0->rgb()[0], v0->rgb()[1], v0->rgb()[2]);
					glBegin(GL_QUADS);
					switch (m_normal_mode)
					{
					case NORMAL_MODE::VERTEX:
						n = v0->normal();
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(v0->point()[0], v0->point()[1], v0->point()[2]);

						n = v1->normal();
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(v1->point()[0], v1->point()[1], v1->point()[2]);

						ratio = (e0->bd_cross() - v2->point()).norm() / (v1->point() - v2->point()).norm();
						n = v2->normal() * ratio + v1->normal() * (1 - ratio);
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(e0->bd_cross()[0], e0->bd_cross()[1], e0->bd_cross()[2]);

						ratio = (e1->bd_cross() - v2->point()).norm() / (v0->point() - v2->point()).norm();
						n = v2->normal() * ratio + v0->normal() * (1 - ratio);
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(e1->bd_cross()[0], e1->bd_cross()[1], e1->bd_cross()[2]);
						break;
					case NORMAL_MODE::FACE:
						n = pF->normal();
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(v0->point()[0], v0->point()[1], v0->point()[2]);
						glVertex3d(v1->point()[0], v1->point()[1], v1->point()[2]);
						glVertex3d(e0->bd_cross()[0], e0->bd_cross()[1], e0->bd_cross()[2]);
						glVertex3d(e1->bd_cross()[0], e1->bd_cross()[1], e1->bd_cross()[2]);
						break;
					};
					glEnd();

					break;

				case 1:
					e0 = bdCrosses[0];
					v0 = m_pMesh->eOppoVertex(e0, pF);
					v1 = m_pMesh->edgeVertex1(e0);
					v2 = m_pMesh->edgeVertex2(e0);
					e1 = m_pMesh->vertexEdge(v0, v2);
					e2 = m_pMesh->vertexEdge(v0, v1);
					oriented = m_pMesh->orientEdges(e1, e2);
					assert(oriented);
					v1 = m_pMesh->eOppoVertex(e1, pF);
					v2 = m_pMesh->eOppoVertex(e2, pF);

					glColor3f(v2->rgb()[0], v2->rgb()[1], v2->rgb()[2]);
					glBegin(GL_TRIANGLES);
					switch (m_normal_mode)
					{
					case NORMAL_MODE::VERTEX:
						n = v2->normal();
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(v2->point()[0], v2->point()[1], v2->point()[2]);

						n = v0->normal();
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(v0->point()[0], v0->point()[1], v0->point()[2]);

						ratio = (e0->bd_cross() - v2->point()).norm() / (v1->point() - v2->point()).norm();
						n = v2->normal() * ratio + v1->normal() * (1 - ratio);
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(e0->bd_cross()[0], e0->bd_cross()[1], e0->bd_cross()[2]);
						break;
					case NORMAL_MODE::FACE:
						n = pF->normal();
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(v2->point()[0], v2->point()[1], v2->point()[2]);
						glVertex3d(v0->point()[0], v0->point()[1], v0->point()[2]);
						glVertex3d(e0->bd_cross()[0], e0->bd_cross()[1], e0->bd_cross()[2]);
						break;
					};
					glEnd();

					glColor3f(v1->rgb()[0], v1->rgb()[1], v1->rgb()[2]);
					glBegin(GL_TRIANGLES);
					switch (m_normal_mode)
					{
					case NORMAL_MODE::VERTEX:
						n = v1->normal();
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(v1->point()[0], v1->point()[1], v1->point()[2]);

						ratio = (e0->bd_cross() - v2->point()).norm() / (v1->point() - v2->point()).norm();
						n = v2->normal() * ratio + v1->normal() * (1 - ratio);
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(e0->bd_cross()[0], e0->bd_cross()[1], e0->bd_cross()[2]);

						n = v0->normal();
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(v0->point()[0], v0->point()[1], v0->point()[2]);
						break;
					case NORMAL_MODE::FACE:
						n = pF->normal();
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(v1->point()[0], v1->point()[1], v1->point()[2]);
						glVertex3d(e0->bd_cross()[0], e0->bd_cross()[1], e0->bd_cross()[2]);
						glVertex3d(v0->point()[0], v0->point()[1], v0->point()[2]);
						break;
					};
					glEnd();

					break;

				case 0:
					int vid = 0;
					int start = 0;
					for (CVDMesh::FaceVertexIterator viter(pF); !viter.end(); viter++)
					{
						v[vid] = *viter;
						if (v[vid]->isomark() == -1)
						{
							start = vid;
						}
						vid++;
					}
					v0 = v[start];
					v1 = v[(start + 1) % 3];
					v2 = v[(start + 2) % 3];

					glColor3f(v0->rgb()[0], v0->rgb()[1], v0->rgb()[2]);
					glBegin(GL_TRIANGLES);
					switch (m_normal_mode)
					{
					case NORMAL_MODE::VERTEX:
						n = v0->normal();
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(v0->point()[0], v0->point()[1], v0->point()[2]);

						n = v1->normal();
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(v1->point()[0], v1->point()[1], v1->point()[2]);

						n = v2->normal();
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(v2->point()[0], v2->point()[1], v2->point()[2]);

						break;
					case NORMAL_MODE::FACE:
						n = pF->normal();
						glNormal3d(a * n[0], a * n[1], a * n[2]);
						glVertex3d(v0->point()[0], v0->point()[1], v0->point()[2]);
						glVertex3d(v1->point()[0], v1->point()[1], v1->point()[2]);
						glVertex3d(v2->point()[0], v2->point()[1], v2->point()[2]);
						break;
					};
					glEnd();
					break;
				}
			}
		}
	}
	
	glBegin(GL_TRIANGLES);

	for (CVDMesh::MeshFaceIterator fiter(m_pMesh); !fiter.end(); ++fiter)
	{
		CViewerFace * pf = *fiter;
		if (pf->boundary_mark() && m_showCut)
		{
			continue;
		}

		if (m_show_face_rgb)
		{
			glColor3d(pf->rgb()[0], pf->rgb()[1], pf->rgb()[2]);
		}
		for (CVDMesh::FaceVertexIterator fviter(pf); !fviter.end(); ++fviter)
		{
			CViewerVertex * v = *fviter;
			CPoint pt = v->point();
			CPoint n;

			switch (m_display_mode)
			{
			case DRAW_MODE::FLAT:
			case DRAW_MODE::FLATLINES:
				m_normal_mode = NORMAL_MODE::FACE;
				break;
			case DRAW_MODE::SMOOTH:
				m_normal_mode = NORMAL_MODE::VERTEX;
				break;
			}

			switch (m_normal_mode)
			{
			case NORMAL_MODE::VERTEX:
				n = v->normal();
				break;
			case NORMAL_MODE::FACE:
				n = pf->normal();
				break;
			};

			/*
			if (show_normal_map)
			{
			n = v->normal_map();
			}
			*/
			if (m_reverse_normal)
				glNormal3d(-n[0], -n[1], -n[2]);
			else
				glNormal3d(n[0], n[1], n[2]);

			CPoint2 uv = v->uv();
			switch (m_texture_mode)
			{
			case TEXTURE_MODE::REPLACE:
			case TEXTURE_MODE::MODULATE:
				glTexCoord2d(uv[0], uv[1]);
				break;
			case TEXTURE_MODE::NONE:
				break;
			}

			// if face colors not specified use vertex colors
			if (!m_show_face_rgb)
			{
				if (m_show_vertex_rgb)
				{
					glColor3d(v->rgb()[0], v->rgb()[1], v->rgb()[2]);
				}
				else
					glColor3d(m_default_rgb[0], m_default_rgb[1], m_default_rgb[2]);
			}

			switch (m_geometry_mode)
			{
			case GEOMETRY_MODE::GEOMETRY:
				glVertex3d(pt[0], pt[1], pt[2]);
				break;
			case GEOMETRY_MODE::PARAMETER:
				glVertex3d(uv[0], uv[1], 0);
				break;
			}
		}
	}
	glEnd();

	m_show_edge = m_display_mode == DRAW_MODE::FLATLINES;

	if (m_show_edge)
	{
		drawEdges();
	}

	if (m_show_boundary)
	{
		drawBoundary(m_pMesh);
	}

	if (m_show_backface)
	{
		drawBackface(m_pMesh);
	}

	if (m_show_sharp_edge)
	{
		drawSharpEdge(m_pMesh);
	}

	if (m_show_picked)
	{
		drawSelectedVertcies();
	}
}

void CMeshViewer::drawBackface(CVDMesh * pMesh)
{
	// Draw Back faces;
	glFrontFace(GL_CW);
	//glCullFace(GL_FRONT);
	//glColor3f(0, 1, 0);
	glColor3f(m_default_rgb[0], m_default_rgb[1], m_default_rgb[2]);

	switch (m_texture_mode)
	{
	case TEXTURE_MODE::REPLACE:
	case TEXTURE_MODE::MODULATE:
		glBindTexture(GL_TEXTURE_2D, m_texName);
		break;
	case TEXTURE_MODE::NONE:
		break;
	}


	glBegin(GL_TRIANGLES);

	for (CVDMesh::MeshFaceIterator fiter(pMesh); !fiter.end(); fiter++)
	{
		CViewerFace * pf = *fiter;

		for (CVDMesh::FaceVertexIterator fviter(pf); !fviter.end(); ++fviter)
		{
			CViewerVertex * v = *fviter;
			CPoint pt = v->point();
			CPoint n;
			switch (m_normal_mode)
			{
			case NORMAL_MODE::FACE:
				n = pf->normal();
				break;
			case NORMAL_MODE::VERTEX:
				n = v->normal();
				break;
			}

			CPoint2 uv = v->uv();
			//glNormal3d( -n[0], -n[1], -n[2] );
			if (m_reverse_normal)
				glNormal3d(n[0], n[1], n[2]);
			else
				glNormal3d(-n[0], -n[1], -n[2]);

			glTexCoord2d(uv[0], uv[1]);

			switch (m_texture_mode)
			{
			case TEXTURE_MODE::REPLACE:
			case TEXTURE_MODE::MODULATE:
				glTexCoord2d(uv[0], uv[1]);
				break;
			case TEXTURE_MODE::NONE:
				break;
			}


			switch (m_geometry_mode)
			{
			case GEOMETRY_MODE::GEOMETRY:
				glVertex3d(pt[0], pt[1], pt[2]);
				break;
			case GEOMETRY_MODE::PARAMETER:
				glVertex3d(uv[0], uv[1], 0);
				break;
			}
		}
	}
	glEnd();
	glFrontFace(GL_CCW);

}



void CMeshViewer::drawWireframe()
{
	if (!m_pMesh) return;

	glPolygonMode(GL_FRONT, GL_LINE);

	//glBindTexture(GL_TEXTURE_2D, texName);
	glBegin(GL_TRIANGLES);

	for (CVDMesh::MeshFaceIterator fiter(m_pMesh); !fiter.end(); ++fiter)
	{
		CViewerFace * pf = *fiter;

		for (CVDMesh::FaceVertexIterator fviter(pf); !fviter.end(); ++fviter)
		{
			CViewerVertex * v = *fviter;
			CPoint pt = v->point();
			CPoint n = v->normal();
			glNormal3d(n[0], n[1], n[2]);
			glColor3d(m_default_rgb[0], m_default_rgb[1], m_default_rgb[2]);

			/*
			if (!show_uv)
			glVertex3d(pt[0], pt[1], pt[2]);
			else
			glVertex3d(uv[0], uv[1], 0);
			*/
			glVertex3d(pt[0], pt[1], pt[2]);
		}
	}
	glEnd();

	glPolygonMode(GL_FRONT, GL_FILL);

	if (m_show_picked)
	{
		drawSelectedVertcies();
	}
}

void CMeshViewer::paintGL()
{
	//std::cout << "paintGL" << std::endl;

	/* clear frame buffer */
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glViewport(0, 0, width(), height());
	// draw gradient background
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(-1, 1, -1, 1);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glDepthMask(GL_FALSE);
	glDisable(GL_TEXTURE_2D);

	glPolygonMode(GL_FRONT, GL_FILL);
	glBegin(GL_QUADS);
	// bottom
	glColor3f(0.5, 0.5, 1.0);
	//glColor3f(1, 1, 1);
	glVertex2f(-1, -1);
	glVertex2f(1, -1);
	// top
	glColor3f(0.0, 0.0, 0.0);
	//glColor3f(1, 1, 1);
	glVertex2f(1, 1);
	glVertex2f(-1, 1);
	glEnd();

	// draw mesh on the background
	glClear(GL_DEPTH_BUFFER_BIT);
	glDepthMask(GL_TRUE);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fovy(), (double)width() / height(), zNear(), zFar());
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//setup GL states

	GLfloat lightOneColor[] = { 1, 1, 1, 1 };
	GLfloat globalAmb[] = { .1f, .1f, .1f, 1.0f };
	GLfloat lightOnePosition[] = { .0f, .0f, 1.0f, 0.0f };

	glEnable(GL_CULL_FACE);
	glFrontFace(GL_CCW);
	//glClearColor(0.1, 0.1, 0.25, 0.5);
	glShadeModel(GL_SMOOTH);

	glEnable(GL_LIGHT1);
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);

	glLightfv(GL_LIGHT1, GL_DIFFUSE, lightOneColor);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, globalAmb);

	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glLightfv(GL_LIGHT1, GL_POSITION, lightOnePosition);

	const GLfloat specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };

	glLightfv(GL_LIGHT1, GL_SPECULAR, specular);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 64.0f);

	GLfloat mat_ambient[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	GLfloat mat_diffuse[] = { 0.01f, 0.01f, 0.01f, 1.0f };
	GLfloat mat_specular[] = { 0.5f, 0.5f, 0.5f, 1.0f };
	GLfloat mat_shininess[] = { 32 };

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

	setupLight();
	/* transform from the eye coordinate system to the world system */
	setupEye();
	glPushMatrix();
	/* transform from the world to the ojbect coordinate system */
	setupObject();

	if (m_polygon_mode == POLYGON_MODE::LINE)
	{
		glPolygonMode(GL_FRONT, GL_LINE);
	}
	switch (m_display_mode)
	{
	case DRAW_MODE::POINTS:
		drawVertices();
		break;
	case DRAW_MODE::WIREFRAME:
		drawWireframe();
		break;
	case DRAW_MODE::TEXTURE:
		if (m_texture_mode != TEXTURE_MODE::NONE)
		{
			glEnable(GL_TEXTURE_2D);
			drawTexture();
		}
		else
		{
			drawMesh();
		}
		break;
	default:
		drawMesh();
	}

	// draw foreground
	if (!m_strokes.curr.empty())
	{
		glClear(GL_DEPTH_BUFFER_BIT);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluOrtho2D(-1, 1, -1, 1);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		glDisable(GL_DEPTH_TEST);
		glDisable(GL_LIGHTING);
		glDepthMask(GL_FALSE);
		glDisable(GL_TEXTURE_2D);

		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

		glLineWidth(5);
		glColor3f(1, 0, 0);
		glBegin(GL_LINES);
		if (!m_strokes.empty())
		{
			for (int i = 0; i < m_strokes.all.size(); i++)
			{
				for (int j = 0; j < m_strokes.all[i].size() - 1; j++)
				{
					double x1 = (double)m_strokes.all[i][j].rx();
					double y1 = height() - (double)m_strokes.all[i][j].ry();
					x1 = x1 * 2 / width() - 1;
					y1 = y1 * 2 / height() - 1;
					double x2 = (double)m_strokes.all[i][j + 1].rx();
					double y2 = height() - (double)m_strokes.all[i][j + 1].ry();
					x2 = x2 * 2 / width() - 1;
					y2 = y2 * 2 / height() - 1;
					glVertex2f(x1, y1);
					glVertex2f(x2, y2);
				}
			}
		}
		for (int j = 0; j < m_strokes.curr.size() - 1; j++)
		{
			double x1 = (double)m_strokes.curr[j].rx();
			double y1 = height() - (double)m_strokes.curr[j].ry();
			x1 = x1 * 2 / width() - 1;
			y1 = y1 * 2 / height() - 1;
			double x2 = (double)m_strokes.curr[j + 1].rx();
			double y2 = height() - (double)m_strokes.curr[j + 1].ry();
			x2 = x2 * 2 / width() - 1;
			y2 = y2 * 2 / height() - 1;
			glVertex2f(x1, y1);
			glVertex2f(x2, y2);
		}
		glEnd();
		glLineWidth(1);

		glDisable(GL_LINE_SMOOTH);
		glBlendFunc(GL_NONE, GL_NONE);
		glDisable(GL_BLEND);
		glEnable(GL_LIGHTING);
	}

	glClear(GL_DEPTH_BUFFER_BIT);
	glDepthMask(GL_TRUE);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fovy(), (double)width() / height(), zNear(), zFar());
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glPopMatrix();
}

void CMeshViewer::resizeGL(int width, int height)
{
	glViewport(0, 0, width, height);
	updateProjectionMatrix();
	updateGL();
}

void CMeshViewer::updateProjectionMatrix()
{
	makeCurrent();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fovy(), (GLdouble)width() / (GLdouble)height(), zNear(), zFar());
}

void CMeshViewer::initializeGL()
{
	//initialize display

	glViewport(0, 0, width(), height());
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fovy(), (double)width() / height(), zNear(), zFar());
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

/*! setup the object, transform from the world to the object coordinate system */
void CMeshViewer::setupObject(void)
{
	double rot[16];

	glTranslated(m_ObjTrans[0], m_ObjTrans[1], m_ObjTrans[2]);
	m_ObjRot.convert(rot);
	glMultMatrixd((GLdouble *)rot);
}

/*! the eye is always fixed at world z = +5 */

void CMeshViewer::setupEye()
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	CQrot rot = m_EyeRot;
	rot = rot ^ (-1);

	double mrot[16];
	m_EyeRot.convert(mrot);
	glMultMatrixd((GLdouble *)mrot);

	//gluLookAt(0, 0, 5, 0, 0, 0, 0, 1, 0);

	glTranslated(-m_EyeTrans[0], -m_EyeTrans[1], -m_EyeTrans[2]);

}

/*! setup light */
void CMeshViewer::setupLight()
{
	CPoint position(0, 0, 1);
	GLfloat lightOnePosition[4] = { position[0], position[1], position[2], 0 };
	glLightfv(GL_LIGHT1, GL_POSITION, lightOnePosition);
	//glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
}

void CMeshViewer::loadMeshFile(const char * meshfile, std::string fileExt)
{
	m_showCut = false;
	m_selection = -1;
	m_part_count = 1;
	if (m_isoline_finder)
	{
		delete m_isoline_finder;
		m_isoline_finder = nullptr;
		m_isoline_sequence.clear();
		m_isoline_metrics.clear();
		m_isoline.clear();
		m_boundary_faces.clear();
	}
	if (m_pMesh)
	{
		delete m_pMesh;
		m_pMesh = nullptr;
		m_picked_vertices.clear();
		m_picked_faces.clear();
	}
	if (!m_strokes.empty())
	{
		m_strokes.clear();
	}
	m_pMesh = new CVDMesh;
	m_isoline_finder = new CIsoline;

	if (fileExt == "obj")
	{
		m_pMesh->read_obj(meshfile);
	}

	if (fileExt == "m")
	{
		m_pMesh->read_m(meshfile);
	}

	if (fileExt == "vef")
	{
		m_pMesh->read_vef(meshfile);
	}

	COperator<CVDMesh> pS(m_pMesh);
	pS._normalize();
	pS._calculate_face_vertex_normal();

	bool has_rgb = false;
	for (CVDMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++)
	{
		CViewerVertex * pV = *viter;
		CPoint dp = pV->rgb() - CPoint(1, 1, 1);
		if (dp.norm() > 0) has_rgb = true;
	}
	if (!has_rgb)
	{
		for (CVDMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++)
		{
			CViewerVertex * pV = *viter;
			pV->rgb() = m_default_rgb;
		}
	}

	std::cout << '\n' << "Mesh imported:" << std::endl;
	std::cout << "Vertices: " << m_pMesh->numVertices();
	std::cout << " Faces: " << m_pMesh->numFaces() << std::endl << std::endl;

	m_isoline_finder->computeMeshInfo(m_pMesh);
	if (m_doNormalFilter)
	{
		normalFilter();
	}
	m_show_picked = true;
}

void CMeshViewer::loadTextureFile(const char * texturefile, std::string fileExt)
{
	if (m_pTexture) delete m_pTexture;

	m_pTexture = new RgbImage;

	if (fileExt == "bmp")
	{
		m_pTexture->LoadBmpFile(texturefile);
		initializeTexture();
	}
}

/*! initialize bitmap image texture */
void CMeshViewer::initializeTexture()
{
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glGenTextures(1, &m_texName);
	glBindTexture(GL_TEXTURE_2D, m_texName);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	int ImageWidth = m_pTexture->GetNumCols();
	int ImageHeight = m_pTexture->GetNumRows();
	GLubyte * ptr = (GLubyte *)m_pTexture->ImageData();

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,
		ImageWidth,
		ImageHeight,
		0,
		GL_RGB,
		GL_UNSIGNED_BYTE,
		ptr);

	switch (m_texture_mode)
	{
	case TEXTURE_MODE::REPLACE:
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		break;
	case TEXTURE_MODE::MODULATE:
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
		break;
	}
}

void CMeshViewer::openMesh()
{
	m_filename = QFileDialog::getOpenFileName(this,
		tr("Open Surface meshes"),
		tr("../models/"),
		tr("obj Files (*.obj);;"
		"m Files (*.m);;"));

	QFileInfo * fileInfo = new QFileInfo(m_filename);
	QString fileExt = fileInfo->suffix();
	std::string sFileExt = fileExt.toStdString();
	if (!m_filename.isEmpty())
	{
		QByteArray byteArray = m_filename.toUtf8();
		const char * _filename = byteArray.constData();
		std::string sFilename = m_filename.toStdString();
		loadMeshFile(_filename, sFileExt);
	}

	updateGL();
}

void CMeshViewer::openTexture()
{
	QString m_filename = QFileDialog::getOpenFileName(this,
		tr("Open Texture images"),
		tr("../textures/"),
		tr("Bmp Files (*.bmp);;"));

	QFileInfo * fileInfo = new QFileInfo(m_filename);
	QString fileExt = fileInfo->suffix();
	std::string sFileExt = fileExt.toStdString();
	if (!m_filename.isEmpty())
	{
		QByteArray byteArray = m_filename.toUtf8();
		const char * _filename = byteArray.constData();
		std::string sFilename = m_filename.toStdString();
		loadTextureFile(_filename, sFileExt);
	}

	updateGL();
}


void CMeshViewer::mousePressEvent(QMouseEvent * mouseEvent)
{
	/* set up an arcball around the Eye's center
	switch y coordinates to right handed system  */


	latestMousePos = mouseEvent->pos();
	gButton = mouseEvent->button();

	if (mouseEvent->button() != Qt::MidButton)
	{
		if (m_isBrushMode)
		{
			if (mouseEvent->button() == Qt::LeftButton)
			{
				m_strokes.curr.clear();
				m_strokes.curr.push_back(latestMousePos);
				pickRay(latestMousePos.rx(), latestMousePos.ry()); //generate the picking ray

				switch (m_brush_mode)
				{
				case BRUSH_MODE::PART:
					if (m_cut_mesh_status == 1)
					{
						for (int id = 0; id < m_isoline_sequence.size(); id++)
						{
							m_isoline_finder->clearMarks(m_isoline_sequence[id]);
						}
						if (!m_strokes.empty())
						{
							m_strokes.part_pop_back();
						}
					}
					m_strokes.curr_ends.start = nullptr;
					m_strokes.curr_ends.end = nullptr;
					/*for (int i = 0; i < m_strokes.curr_faces.size(); i++)
					{
					CViewerFace * pF = m_strokes.curr_faces[i];
					pF->picked() = false;
					}*/
					m_strokes.curr_faces.clear();
					m_isoline_sequence.clear();

					if (!m_updateMode)
					{
						m_strokes.part_clear();
						for (int i = 0; i < m_picked_faces.size(); i++)
						{
							CViewerFace * pF = m_picked_faces[i];
							pF->picked() = false;
						}
						m_picked_faces.clear();
					}

					if (m_pMesh)
					{
						CViewerVertex * pV = pickVertex(); //return the selected vertex
						if (pV)
						{
							m_strokes.curr_ends.start = pV;
						}
						CViewerFace * pF = pickFace();
						if (pF)
						{
							pF->picked() = true;
							m_strokes.curr_faces.push_back(pF);
						}
					}
					break;
				case BRUSH_MODE::PATCH:
					m_strokes.curr_end_faces.start = nullptr;
					m_strokes.curr_end_faces.end = nullptr;

					if (m_pMesh)
					{
						CViewerFace * pF = pickFace();
						if (pF)
						{
							pF->picked() = true;
							m_strokes.curr_end_faces.start = pF;
						}
					}
					break;
				default:
					break;
				}
				
				makeCurrent();
				updateGL();
			}
		}
		if (m_isSelectionMode && m_pMesh)
		{
			pickRay(latestMousePos.rx(), latestMousePos.ry()); //generate the picking ray 
			CViewerVertex *pV = nullptr;
			CViewerFace   *pF = nullptr;

			switch (m_pick_mode)
			{
			case PICK_MODE::VERTEX:

				pV = pickVertex(); //return the selected vertex
				if (pV)
				{
					if (mouseEvent->button() == Qt::LeftButton)
					{
						if (!pV->picked())
						{
							pV->picked() = true;
							m_picked_vertices.push_back(pV);
						}
					}
					if (mouseEvent->button() == Qt::RightButton && pV->picked())
					{
						pV->picked() = false;
						VL::iterator iter, iter0, iter1;
						iter0 = m_picked_vertices.begin();
						iter1 = m_picked_vertices.end();
						iter = std::find(iter0, iter1, pV);
						m_picked_vertices.erase(iter);
					}
				}
				break;
			case PICK_MODE::FACE:
				pF = pickFace();
				if (pF)
				{
					if (mouseEvent->button() == Qt::LeftButton)
					{
						if (!pF->picked())
						{
							pF->picked() = true;
							m_picked_faces.push_back(pF);
						}
					}
					if (mouseEvent->button() == Qt::RightButton)
					{
						pF->picked() = false;
						FL::iterator iter, iter0, iter1;
						iter0 = m_picked_faces.begin();
						iter1 = m_picked_faces.end();
						iter = std::find(iter0, iter1, pF);
						m_picked_faces.erase(iter);
					}
				}
				break;
			}
			makeCurrent();
			updateGL();
		}
		arcball = CArcball(width(), height(), latestMousePos.rx() - width() / 2, height() - latestMousePos.ry() - height() / 2);
	}
}

void CMeshViewer::mouseMoveEvent(QMouseEvent * mouseEvent)
{
	CPoint trans;
	CQrot       rot;

	makeCurrent();

	/* rotation, call arcball */
	if (gButton == Qt::LeftButton)
	{
		latestMousePos = mouseEvent->pos();
		if (m_isBrushMode)
		{
			m_strokes.curr.push_back(latestMousePos);
			arcball = CArcball(width(), height(), latestMousePos.rx() - width() / 2, height() - latestMousePos.ry() - height() / 2);
		}
		else if (m_isSelectionMode && m_pMesh)
		{
			pickRay(latestMousePos.rx(), latestMousePos.ry()); //generate the picking ray 
			CViewerVertex *pV = nullptr;
			CViewerFace   *pF = nullptr;

			switch (m_pick_mode)
			{
			case PICK_MODE::VERTEX:

				pV = pickVertex(); //return the selected vertex
				if (pV)
				{
					if (!pV->picked())
					{
						pV->picked() = true;
						m_picked_vertices.push_back(pV);
					}
				}
				break;
			case PICK_MODE::FACE:
				pF = pickFace();
				if (pF && !pF->picked())
				{
					pF->picked() = true;
					m_picked_faces.push_back(pF);
				}
				break;
			}
			//makeCurrent();
			//updateGL();
			arcball = CArcball(width(), height(), latestMousePos.rx() - width() / 2, height() - latestMousePos.ry() - height() / 2);
		}
		else
		{
			rot = arcball.update(mouseEvent->pos().rx() - width() / 2, height() - mouseEvent->pos().ry() - height() / 2);
			switch (m_manipulation_mode)
			{
			case MANIPULATION_MODE::EYE:
				m_EyeRot = m_EyeRot * rot;
				break;
			case MANIPULATION_MODE::OBJECT:
				m_ObjRot = rot * m_ObjRot;
				break;
			}
		}

	}

	/*xy translation */
	if (gButton == Qt::MidButton)
	{
		double scale = 10. / height();
		trans = CPoint(scale*(mouseEvent->pos().rx() - latestMousePos.rx()), scale*(latestMousePos.ry() - mouseEvent->pos().ry()), 0);
		latestMousePos = mouseEvent->pos();

		switch (m_manipulation_mode)
		{
		case MANIPULATION_MODE::EYE:
			m_EyeTrans = m_EyeTrans - trans;
			break;
		case MANIPULATION_MODE::OBJECT:
			m_ObjTrans = m_ObjTrans + trans;
			break;
		}
	}

	if (gButton == Qt::RightButton)
	{
		if (m_isSelectionMode && m_pMesh)
		{
			latestMousePos = mouseEvent->pos();
			pickRay(latestMousePos.rx(), latestMousePos.ry()); //generate the picking ray 
			CViewerVertex *pV = nullptr;
			CViewerFace   *pF = nullptr;

			switch (m_pick_mode)
			{
			case PICK_MODE::VERTEX:

				pV = pickVertex(); //return the selected vertex
				if (pV && pV->picked())
				{
					pV->picked() = false;
					VL::iterator iter, iter0, iter1;
					iter0 = m_picked_vertices.begin();
					iter1 = m_picked_vertices.end();
					iter = std::find(iter0, iter1, pV);
					m_picked_vertices.erase(iter);
				}
				break;
			case PICK_MODE::FACE:
				pF = pickFace();
				if (pF)
				{
					pF->picked() = false;
					FL::iterator iter, iter0, iter1;
					iter0 = m_picked_faces.begin();
					iter1 = m_picked_faces.end();
					iter = std::find(iter0, iter1, pF);
					m_picked_faces.erase(iter);
				}
				break;
			}
			//makeCurrent();
			//updateGL();
			arcball = CArcball(width(), height(), latestMousePos.rx() - width() / 2, height() - latestMousePos.ry() - height() / 2);
		}
		else
			/* zoom in and out */
		{
			double scale = 10. / height();
			trans = CPoint(0, 0, scale*(latestMousePos.ry() - mouseEvent->pos().ry()));
			latestMousePos = mouseEvent->pos();
			switch (m_manipulation_mode)
			{
			case MANIPULATION_MODE::EYE:
				m_EyeTrans = m_EyeTrans - trans;
				break;
			case MANIPULATION_MODE::OBJECT:
				m_ObjTrans = m_ObjTrans + trans;
				break;
			}
		}
	}

	// update OpenGL, trigger re-draw
	updateGL();
}

void CMeshViewer::mouseReleaseEvent(QMouseEvent * mouseEvent)
{
	latestMousePos = mouseEvent->pos();
	if (mouseEvent->button() == Qt::LeftButton)
	{
		if (m_isBrushMode)
		{
			//m_strokes.curr.push_back(latestMousePos);
			pickRay(latestMousePos.rx(), latestMousePos.ry()); //generate the picking ray
			switch (m_brush_mode)
			{
			case BRUSH_MODE::PART:
				if (m_pMesh)
				{
					CViewerVertex * pV = pickVertex(); //return the selected vertex
					if (pV)
					{
						m_strokes.curr_ends.end = pV;
					}
					CViewerFace * pF = pickFace();
					if (pF && !pF->picked())
					{
						pF->picked() = true;
						m_strokes.curr_faces.push_back(pF);
					}
					if (m_strokes.curr_faces.size() == 2 && m_strokes.curr_ends.start && m_strokes.curr_ends.end)
					{
						m_strokes.curr_faces = findPath(m_strokes.curr_faces.front(), m_strokes.curr_faces.back());
						m_strokes.part_push_back();
						if (m_updateMode && m_part_count > 1)
						{
							for (int i = 0; i < m_boundary_faces.back().size(); i++)
							{
								CViewerFace *pF = m_boundary_faces.back()[i];
								pF->boundary_mark() = false;
								for (CVDMesh::FaceHalfedgeIterator hiter(pF); !hiter.end(); hiter++)
								{
									CViewerEdge *pE = m_pMesh->halfedgeEdge(*hiter);
									pE->bd_cross() = CPoint(0, 0, 0);
								}
							}
							m_boundary_faces.pop_back();
							m_part_count--;
						}
						if (!m_slowMotion)
						{
							computeField();
							int failed = computeIsolines();
							if (failed)
							{
								for (int id = 0; id < m_isoline_sequence.size(); id++)
								{
									m_isoline_finder->clearMarks(m_isoline_sequence[id]);
								}
								m_isoline_sequence.clear();
								m_strokes.part_pop_back();
							}
							else
							{
								m_cut_mesh_status = cutMesh();
							}
						}
					}
					else
					{
						QMessageBox box(QMessageBox::Warning, "Insufficient constraints", "Please redraw a stroke to guarantee 2 points picked.", QMessageBox::Ok, NULL);
						box.exec();
					}
				}
				break;
			case BRUSH_MODE::PATCH:
				if (m_pMesh)
				{
					CViewerFace * pF = pickFace();
					if (pF && !pF->picked())
					{
						pF->picked() = true;
						m_strokes.curr_end_faces.end = pF;
						if (m_strokes.curr_end_faces.start)
						{
							m_strokes.patch_push_back();
							if (!m_slowMotion)
							{
								regionGrow();
							}
						}
					}
				}
				break;
			default:
				break;
			}
		}
	}
	gButton = Qt::NoButton;
	updateGL();
}

void CMeshViewer::wheelEvent(QWheelEvent * mouseEvent)
{
	// scroll the wheel to scale the view port
	double moveAmount = (double)mouseEvent->delta() / (120.0*8.0);
	CPoint trans = CPoint(0.0, 0.0, moveAmount);
	switch (m_manipulation_mode)
	{
	case MANIPULATION_MODE::EYE:
		m_EyeTrans = m_EyeTrans - trans;
		break;
	case MANIPULATION_MODE::OBJECT:
		m_ObjTrans = m_ObjTrans + trans;
		break;
	}
	updateGL();
	mouseEvent->accept();
}

void CMeshViewer::drawVertices()
{
	if (!m_pMesh) return;

	glDisable(GL_LIGHTING);
	glPointSize(1);
	glColor3f(1, 1, 1);
	glBegin(GL_POINTS);

	for (CVDMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++)
	{
		CViewerVertex * pV = *viter;

		CPoint p = pV->point();
		CPoint2 uv = pV->uv();

		switch (m_geometry_mode)
		{
		case GEOMETRY_MODE::GEOMETRY:
			glVertex3d(p[0], p[1], p[2]);
			break;
		case GEOMETRY_MODE::PARAMETER:
			glVertex3d(uv[0], uv[1], 0);
			break;
		}
	}
	glEnd();
	glEnable(GL_LIGHTING);
}

void CMeshViewer::drawEdges()
{
	if (!m_pMesh) return;

	glDisable(GL_LIGHTING);
	glColor3f(0, 0, 0);
	glBegin(GL_LINES);

	for (CVDMesh::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter)
	{
		CViewerEdge   * pE = *eiter;

		if (m_show_sharp_edge && pE->sharp())
			continue;

		CViewerVertex * v1 = m_pMesh->edgeVertex1(pE);
		CViewerVertex * v2 = m_pMesh->edgeVertex2(pE);

		CPoint p1 = v1->point();
		CPoint p2 = v2->point();

		CPoint2 uv1 = v1->uv();
		CPoint2 uv2 = v2->uv();

		switch (m_geometry_mode)
		{
		case GEOMETRY_MODE::GEOMETRY:
			glVertex3d(p1[0], p1[1], p1[2]);
			glVertex3d(p2[0], p2[1], p2[2]);
			break;
		case GEOMETRY_MODE::PARAMETER:
			glVertex3d(uv1[0], uv1[1], 0);
			glVertex3d(uv2[0], uv2[1], 0);
			break;
		}
	}
	glEnd();
	glEnable(GL_LIGHTING);
}

void CMeshViewer::showPoints()
{
	m_display_mode = DRAW_MODE::POINTS;
	glDisable(GL_TEXTURE_2D);
	updateGL();
}

void CMeshViewer::showSmooth()
{
	m_display_mode = DRAW_MODE::SMOOTH;
	glDisable(GL_TEXTURE_2D);
	updateGL();
}

void CMeshViewer::showFlat()
{
	m_display_mode = DRAW_MODE::FLAT;
	glDisable(GL_TEXTURE_2D);
	updateGL();
}

void CMeshViewer::showFlatlines()
{
	m_display_mode = DRAW_MODE::FLATLINES;
	glDisable(GL_TEXTURE_2D);
	updateGL();
}

void CMeshViewer::showBoundary()
{
	m_show_boundary = !m_show_boundary;
	updateGL();
}

void CMeshViewer::showWireframe()
{
	m_display_mode = DRAW_MODE::WIREFRAME;
	glDisable(GL_TEXTURE_2D);
	updateGL();
}

void CMeshViewer::showTexture()
{
	m_display_mode = DRAW_MODE::TEXTURE;
	m_texture_mode = TEXTURE_MODE::MODULATE;
	glEnable(GL_TEXTURE_2D);
	updateGL();
}

void CMeshViewer::drawTexture()
{
	if (!m_pMesh) return;
	if (!m_pTexture) return;

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0, 1.0);

	glBindTexture(GL_TEXTURE_2D, m_texName);

	glBegin(GL_TRIANGLES);

	for (CVDMesh::MeshFaceIterator fiter(m_pMesh); !fiter.end(); ++fiter)
	{
		CViewerFace * pf = *fiter;

		for (CVDMesh::FaceVertexIterator fviter(pf); !fviter.end(); ++fviter)
		{
			CViewerVertex * v = *fviter;
			CPoint pt = v->point();
			CPoint n = v->normal();
			CPoint2 uv = v->uv();

			glNormal3d(n[0], n[1], n[2]);


			glTexCoord2d(uv[0], uv[1]);
			glColor3d(m_default_rgb[0], m_default_rgb[1], m_default_rgb[2]);

			/*
			if (!show_uv)
			glVertex3d(pt[0], pt[1], pt[2]);
			else
			glVertex3d(uv[0], uv[1], 0);
			*/
			glVertex3d(pt[0], pt[1], pt[2]);
		}
	}
	glEnd();
}

void CMeshViewer::set_texture_mode(TEXTURE_MODE mode)
{
	m_texture_mode = mode;
	switch (mode)
	{
	case TEXTURE_MODE::REPLACE:
		glEnable(GL_TEXTURE_2D);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		break;
	case TEXTURE_MODE::MODULATE:
		glEnable(GL_TEXTURE_2D);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
		break;
	case TEXTURE_MODE::NONE:
		glDisable(GL_TEXTURE_2D);
		break;
	}
	updateGL();
	makeCurrent();
}

void CMeshViewer::set_polygon_mode(POLYGON_MODE mode)
{
	m_polygon_mode = mode;

	switch (mode)
	{
	case POLYGON_MODE::FILL:
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		break;
	case POLYGON_MODE::LINE:
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		break;
	}
	updateGL();
	makeCurrent();
}

void CMeshViewer::set_orientation_mode(ORIENTATION_MODE mode)
{
	switch (mode)
	{
	case ORIENTATION_MODE::CCW:
		glFrontFace(GL_CCW);
		break;
	case ORIENTATION_MODE::CLW:
		glFrontFace(GL_CW);
		break;
	}
	updateGL();
	makeCurrent();
}

void CMeshViewer::set_edge_flag(bool flag)
{
	m_show_edge = flag;
	updateGL();
	makeCurrent();
}

void CMeshViewer::set_boundary_flag(bool flag)
{
	m_show_boundary = flag;
	updateGL();
	makeCurrent();
}

void CMeshViewer::set_sharp_edge_flag(bool flag)
{
	m_show_sharp_edge = flag;
	updateGL();
	makeCurrent();
}

void CMeshViewer::set_vertex_rgb_flag(bool flag)
{
	m_show_vertex_rgb = flag;
	updateGL();
	makeCurrent();
}

void CMeshViewer::set_face_rgb_flag(bool flag)
{
	m_show_face_rgb = flag;
	updateGL();
	makeCurrent();
}

void CMeshViewer::set_show_cut_flag(bool flag)
{
	m_showCut = flag;
	updateGL();
	makeCurrent();
}

void CMeshViewer::set_backface_flag(bool flag)
{
	m_show_backface = flag;
	updateGL();
	makeCurrent();
}

void CMeshViewer::set_filter_normal_flag(bool flag)
{
	m_doNormalFilter = flag;
	updateGL();
	makeCurrent();
}

void CMeshViewer::set_reverse_normal_flag(bool flag)
{
	m_reverse_normal = flag;
	updateGL();
	makeCurrent();
}

void CMeshViewer::set_picked_flag(bool flag)
{
	m_show_picked = flag;
	updateGL();
	makeCurrent();
}

void CMeshViewer::set_slow_motion_flag(bool flag)
{
	m_slowMotion = flag;
	updateGL();
	makeCurrent();
}

void CMeshViewer::set_update_mode(bool flag)
{
	m_updateMode = flag;
	if (!m_updateMode)
	{
		m_strokes.part_clear();
		for (int i = 0; i < m_picked_faces.size(); i++)
		{
			CViewerFace * pF = m_picked_faces[i];
			pF->picked() = false;
		}
		m_picked_faces.clear();
	}
}

void CMeshViewer::drawBoundary(CVDMesh * pMesh)
{
	glDisable(GL_LIGHTING);
	glColor3f(0, 0, 1);
	glLineWidth(3.0);
	glBegin(GL_LINES);
	for (CVDMesh::MeshEdgeIterator eiter(pMesh); !eiter.end(); eiter++)
	{
		CViewerEdge * pE = *eiter;
		if (!pE->boundary()) continue;

		CViewerVertex * pV0 = pMesh->edgeVertex1(pE);
		CViewerVertex * pV1 = pMesh->edgeVertex2(pE);

		CPoint p0 = pV0->point();
		CPoint p1 = pV1->point();

		CPoint2 uv0 = pV0->uv();
		CPoint2 uv1 = pV1->uv();

		switch (m_geometry_mode)
		{
		case GEOMETRY_MODE::GEOMETRY:
			glVertex3d(p0[0], p0[1], p0[2]);
			glVertex3d(p1[0], p1[1], p1[2]);
			break;
		case GEOMETRY_MODE::PARAMETER:
			glVertex3d(uv0[0], uv0[1], 0);
			glVertex3d(uv1[0], uv1[1], 0);
			break;
		}

	}
	glEnd();
	glLineWidth(1.0);
	glEnable(GL_LIGHTING);
}

void CMeshViewer::drawSharpEdge(CVDMesh * pMesh)
{
	glDisable(GL_LIGHTING);
	glColor3f(1, 1, 0);
	glLineWidth(3.0);
	glBegin(GL_LINES);
	for (CVDMesh::MeshEdgeIterator eiter(pMesh); !eiter.end(); eiter++)
	{
		CViewerEdge * pE = *eiter;
		if (!pE->sharp()) continue;

		CViewerVertex * pV0 = pMesh->edgeVertex1(pE);
		CViewerVertex * pV1 = pMesh->edgeVertex2(pE);

		CPoint p0 = pV0->point();
		CPoint p1 = pV1->point();

		CPoint2 uv0 = pV0->uv();
		CPoint2 uv1 = pV1->uv();

		switch (m_geometry_mode)
		{
		case GEOMETRY_MODE::GEOMETRY:
			glVertex3d(p0[0], p0[1], p0[2]);
			glVertex3d(p1[0], p1[1], p1[2]);
			break;
		case GEOMETRY_MODE::PARAMETER:
			glVertex3d(uv0[0], uv0[1], 0);
			glVertex3d(uv1[0], uv1[1], 0);
			break;
		}

	}
	glEnd();
	glLineWidth(1.0);
	glEnable(GL_LIGHTING);
}


//copy frame buffer to an image
/*! save frame buffer to an image "snap_k.bmp"
*/
void CMeshViewer::snapShot(std::string file, std::string ext)
{
	static int id = 0;

	int win_width = width();
	int win_height = height();

	GLfloat * buffer = new GLfloat[win_width * win_height * 3];
	assert(buffer);
	glReadBuffer(GL_FRONT_LEFT);
	glReadPixels(0, 0, win_width, win_height, GL_RGB, GL_FLOAT, buffer);

	RgbImage  image(win_height, win_width);

	for (int i = 0; i < win_height; i++)
	for (int j = 0; j < win_width; j++)
	{
		float r = buffer[(i*win_width + j) * 3 + 0];
		float g = buffer[(i*win_width + j) * 3 + 1];
		float b = buffer[(i*win_width + j) * 3 + 2];

		image.SetRgbPixelf(i, j, r, g, b);
	}
	delete[]buffer;

	char name[256];
	std::ostringstream os(name);
	//os << file << "." << ext;
	os << file;
	image.WriteBmpFile(os.str().c_str());

}

void CMeshViewer::saveImage()
{
	QString m_filename = QFileDialog::getSaveFileName(this,
		tr("Save snapshot"),
		tr("../models/"),
		tr("bmp Files (*.bmp);;"
		"m Files (*.m);;"));

	QFileInfo * fileInfo = new QFileInfo(m_filename);
	QString fileExt = fileInfo->suffix();
	std::string sFileExt = fileExt.toStdString();
	if (!m_filename.isEmpty())
	{
		QByteArray byteArray = m_filename.toUtf8();
		//const char * _filename = byteArray.constData();
		std::string sFilename = m_filename.toStdString();
		snapShot(sFilename, sFileExt);
	}

	updateGL();
}


//get the picking ray
void CMeshViewer::pickRay(double x, double y)
{
	GLint viewport[4];
	GLdouble modelview[16], projection[16];
	GLdouble winX, winY, posX, posY, posZ;

	glPushMatrix();
	setupEye();
	setupObject();

	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glGetIntegerv(GL_VIEWPORT, viewport);
	glPopMatrix();

	winX = (double)x;
	winY = height() - (double)y;
	//winY = viewport[3] - (double)y;

	//get the point on the near plane
	gluUnProject(winX, winY, 0.0, modelview, projection, viewport, &posX, &posY, &posZ);
	// std::cout << "1:\t" << posX << " " << posY << " " << posZ << std::endl;
	m_ray[0][0] = posX, m_ray[0][1] = posY, m_ray[0][2] = posZ;

	//get the point on the far plane
	gluUnProject(winX, winY, 1.0, modelview, projection, viewport, &posX, &posY, &posZ);
	// std::cout << "2:\t" << posX << " " << posY << " " << posZ << std::endl;
	m_ray[1][0] = posX, m_ray[1][1] = posY, m_ray[1][2] = posZ;
}

//pick the selected vertex
CViewerVertex * CMeshViewer::pickVertex()
{
	//method one, using the angle of two lines
	double max_cos_theta = DBL_MIN;
	CViewerVertex * v = nullptr;
	QVector3D target = m_ray[1], source = m_ray[0];
	for (CVDMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++)
	{
		CViewerVertex * pV = *viter;
		QVector3D p(pV->point()[0], pV->point()[1], pV->point()[2]);

		double cos_theta = QVector3D::dotProduct(p - source, target - source) / ((p - source).length() * (target - source).length());
		if (cos_theta > EPSILON_COS_THETA_FOR_SELECT_VERTEX && max_cos_theta < cos_theta)
		{
			QVector3D normal(pV->normal()[0], pV->normal()[1], pV->normal()[2]);
			if (QVector3D::dotProduct(normal, target - source) < 0)
			{
				max_cos_theta = cos_theta;
				v = pV;
			}
		}
	}

	//std::cout << "max_cos_theta:" << max_cos_theta << std::endl;
	if (v)
	{
		std::cout << "Picked Vertex ID : " << v->sid() << std::endl;
		//std::cout << v->isomark() << std::endl;
	}

	return v;
}

//pick the selected face
CViewerFace * CMeshViewer::pickFace()
{
	QVector3D p0 = m_ray[0], p1 = m_ray[1];
	//CPoint source(s.x(), s.y(), s.z()), target(t.x(), t.y(), t.z());
	for (CVDMesh::MeshFaceIterator fiter(m_pMesh); !fiter.end(); fiter++)
	{
		CViewerFace * pF = *fiter;
		CPoint _v0 = pF->halfedge()->source()->point();
		CPoint _v1 = pF->halfedge()->he_next()->source()->point();
		CPoint _v2 = pF->halfedge()->he_prev()->source()->point();
		QVector3D v0(_v0[0], _v0[1], _v0[2]);
		QVector3D v1(_v1[0], _v1[1], _v1[2]);
		QVector3D v2(_v2[0], _v2[1], _v2[2]);

		//0. check whether the normal share the same direction with the ray
		QVector3D normal = QVector3D::crossProduct(v1 - v0, v2 - v1).normalized();
		if (QVector3D::dotProduct(normal, p0 - p1) <= 0)
			continue;

		//1. compute the intersection
		QVector3D u = p1 - p0;
		QVector3D w = p0 - v0;

		double D = QVector3D::dotProduct(normal, u);
		double N = -QVector3D::dotProduct(normal, w);
		if (fabs(D) < EPSILON)
			continue;

		double sI = N / D;
		if (sI < 0.0 || sI > 1)
			continue;
		QVector3D I = p0 + sI * u; // the intersection

		//2. check whether the intersection lies in the triangle.
		double theta0 = acos(QVector3D::dotProduct(v0 - I, v1 - I) / ((v0 - I).length() * (v1 - I).length()));
		double theta1 = acos(QVector3D::dotProduct(v1 - I, v2 - I) / ((v1 - I).length() * (v2 - I).length()));
		double theta2 = acos(QVector3D::dotProduct(v2 - I, v0 - I) / ((v2 - I).length() * (v0 - I).length()));
		double total = theta0 + theta1 + theta2;
		if (fabs(total - 2 * PI) < EPSILON)
			return pF;
	}
	return nullptr;
}

void CMeshViewer::enterSelectionMode()
{
	m_isSelectionMode = true;
	m_isBrushMode = false;
	m_strokes.clear();
	updateGL();
	//std::cout << "Enter Selection Mode" << std::endl;
}

void CMeshViewer::quitSelectionMode()
{
	m_isSelectionMode = false;
	//std::cout << "Exit Selection Mode" << std::endl;
}

void CMeshViewer::enterBrushMode()
{
	m_isBrushMode = true;
	m_isSelectionMode = false;
}

void CMeshViewer::quitBrushMode()
{
	m_isBrushMode = false;
	m_strokes.clear();
	updateGL();
}

void CMeshViewer::drawSelectedVertcies()
{

	glDisable(GL_LIGHTING);

	if (!m_isBrushMode)
	{
		glColor3f(0.0, 1.0, 0.0);
		glBegin(GL_TRIANGLES);

		for (int i = 0; i < m_picked_faces.size(); i++)
		{
			CViewerFace * pF = m_picked_faces[i];
			CPoint n = pF->normal();
			for (CVDMesh::FaceVertexIterator fviter(pF); !fviter.end(); fviter++)
			{
				CViewerVertex * pV = *fviter;
				CPoint p = pV->point() + n * 0.001;
				glVertex3d(p[0], p[1], p[2]);
			}
		}
		glEnd();
	}


	glPointSize(10);
	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_NOTEQUAL, 0);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

	glBegin(GL_POINTS);
	for (int i = 0; i < m_picked_vertices.size(); i++)
	{
		CViewerVertex *pV = m_picked_vertices[i];
		CPoint pt = pV->point();
		if (m_geometry_mode == GEOMETRY_MODE::PARAMETER)
		{
			pt = CPoint(pV->uv()[0], pV->uv()[1], 0);
		}

		glColor3f(pV->rgb()[0], pV->rgb()[1], pV->rgb()[2]);
		glVertex3d(pt[0], pt[1], pt[2]);
	}
	glEnd();
	glPointSize(1);
	glDisable(GL_POINT_SMOOTH);

	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	if (!m_isoline_sequence.empty() && m_brush_mode == BRUSH_MODE::PART)
	{
		glColor3f(1, 0, 0);
		glLineWidth(2);
		for (int k = 0; k < m_isoline_sequence.size(); k++)
		{
			glBegin(GL_LINES);
			for (int i = 0; i < m_isoline_sequence[k].size(); i++)
			{
				size_t size_i = m_isoline_sequence[k][i].size();
				for (int j = 0; j < size_i; j++)
				{
					CPoint p1 = m_isoline_sequence[k][i].points[j];
					CPoint p2 = m_isoline_sequence[k][i].points[(j + 1) % size_i];
					glVertex3d(p1[0], p1[1], p1[2]);
					glVertex3d(p2[0], p2[1], p2[2]);
				}
			}
			glEnd();
		}
		glLineWidth(1);
	}

	glDisable(GL_LINE_SMOOTH);
	glBlendFunc(GL_NONE, GL_NONE);
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
}

void CMeshViewer::set_pick_mode(PICK_MODE mode)
{
	m_pick_mode = mode;
	updateGL();
	makeCurrent();
}


void CMeshViewer::set_brush_mode(BRUSH_MODE mode)
{
	m_brush_mode = mode;
	updateGL();
	makeCurrent();
}

void CMeshViewer::loadFromMainWindow(std::string fileName, std::string fExt)
{
	m_filename.push_back(QString::fromStdString(fileName));
	const char * _filename = fileName.c_str();
	if (fExt == "bmp")
		loadTextureFile(_filename, fExt);
	else
		loadMeshFile(_filename, fExt);
}

void CMeshViewer::showField()
{
	if (m_slowMotion)
	{
		computeField();

		// use color to represent scalar field
		double fieldMax = -INFI, fieldMin = INFI, fieldDiam = INFI;
		for (CVDMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++)
		{
			CViewerVertex *pV = *viter;
			if (pV->field() > fieldMax)
			{
				fieldMax = pV->field();
			}
			if (pV->field() < fieldMin)
			{
				fieldMin = pV->field();
			}
		}
		fieldDiam = fieldMax - fieldMin;
		for (CVDMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++)
		{
			CViewerVertex *pV = *viter;
			double f = (pV->field() - fieldMin) / fieldDiam;
			double a = (1 - f) / 0.25;
			int X = floor(a);
			int Y = floor(255 * (a - X));
			int r, g, b;
			switch (X)
			{
			case 0: r = 255; g = Y; b = 0; break;
			case 1: r = 255 - Y; g = 255; b = 0; break;
			case 2: r = 0; g = 255; b = Y; break;
			case 3: r = 0; g = 255 - Y; b = 255; break;
			case 4: r = 0; g = 0; b = 255; break;
			}
			pV->rgb_bak() = pV->rgb();
			pV->rgb() = CPoint(r, g, b) / 255;
		}
		m_show_picked = false;
		updateGL();
	}
}

void CMeshViewer::showIsolines()
{
	if (m_slowMotion)
	{
		int failed = computeIsolines();
		if (failed)
		{
			for (int id = 0; id < m_isoline_sequence.size(); id++)
			{
				m_isoline_finder->clearMarks(m_isoline_sequence[id]);
			}
			m_isoline_sequence.clear();

			m_strokes.part_pop_back();
		}
		m_show_picked = true;
		updateGL();
	}
}

void CMeshViewer::showCutResult()
{
	if (m_slowMotion)
	{
		for (CVDMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++)
		{
			CViewerVertex *pV = *viter;
			pV->rgb() = pV->rgb_bak();
		}
		m_cut_mesh_status = cutMesh();
		updateGL();
	}
}

void CMeshViewer::computeField()
{
	if (!m_strokes.all.empty())
	{
		m_isoline_finder->generateField(m_pMesh, m_strokes.all_ends);
	}
}

int CMeshViewer::computeIsolines()
{
	m_isoline_metrics.clear();
	// N = 15 here
	int N = m_isoline_set_size;
	for (int id = 0; id < N; id++)
	{
		int i = id + 1;
		double thresh = 1.0 * i / (N + 1);
		m_isoline = m_isoline_finder->getIsoline(m_pMesh, thresh);
		filterIsolines();
	}
	if (m_isoline_sequence.size() < N)
	{
		QMessageBox box(QMessageBox::Warning, "Stroke too short", "Please draw a longer stroke.", QMessageBox::Ok, NULL);
		box.exec();
		return 1;
	}
	// compute the metric of each isoline
	ARR r(N, 0), C(N, 0), Delta(N, 0), M(N, 0);
	double t = N / 2.0;
	double sigma = 2;
	double minMetric = 0;
	int selected = 0;
	for (int id = 0; id < N; id++)
	{
		int i = id + 1;
		// centerness
		C[id] = exp(-(i - t) * (i - t) / (2 * t * t));
		// local radius
		r[id] = m_isoline_sequence[id][0].length / (2 * PI);
	}
	for (int id = 0; id < N; id++)
	{
		double sum_f = 0;
		for (int k = 1; k < t; k++)
		{
			double fk = exp(-(k - 1) * (k - 1) / (2 * sigma * sigma));
			int left = id - k >= 0 ? id - k : 0;
			int right = id + k < N ? id + k : N;
			Delta[id] += fk * (2 * r[id] - r[left] - r[right]);
			sum_f += fk;
		}
		// concaveness
		Delta[id] /= sum_f;
	}
	for (int id = 0; id < N; id++)
	{
		M[id] = C[id] * Delta[id];
		std::cout << '#' << id + 1 << '\t' << "r : " << r[id] << '\t' << " Delta : " << Delta[id] << '\t' << "M : " << M[id] << std::endl;
		m_isoline_metrics.push_back(M[id]);
		// select the isoline with lowest score
		if (m_isoline_metrics.back() < minMetric)
		{
			minMetric = m_isoline_metrics.back();
			selected = id;
		}
	}
	m_selection = selected + 1;
	m_isoline.thresh = 1.0 * (selected + 1) / (N + 1);
		
	return 0;
}

void CMeshViewer::filterIsolines()
{
	std::vector<int> strokeIsolineId;
	if (!m_strokes.all.empty() && !m_isoline.empty())
	{
		// select the isolines passing all strokes and not passing existing boundaries
		bool valid = false, inflag = false;
		for (int i = 0; i < m_isoline.size(); i++)
		{
			valid = true;
			for (int j = 0; j < m_strokes.all_faces.size(); j++)
			{
				inflag = false;
				for (int k = 0; k < m_strokes.all_faces[j].size(); k++)
				{
					CViewerFace *pF = m_strokes.all_faces[j][k];
					if (pF->isomark() == i)
					{
						inflag = true;
						break;
					}
				}
				if (!inflag)
				{
					// isoline i doesn't pass stroke j
					valid = false;
					break;
				}
			}
			if (valid)
			{
				strokeIsolineId.push_back(i);
			}
		}

		if (!strokeIsolineId.empty())
		{
			ISO validIsolines;
			for (int i = 0; i < strokeIsolineId.size(); i++)
			{
				validIsolines.push_back(m_isoline[strokeIsolineId[i]]);
			}
			// exclude isolines not passing the stroke
			for (int i = 0; i < m_isoline.size(); i++)
			{
				valid = false;
				for (int j = 0; j < strokeIsolineId.size(); j++)
				{
					if (strokeIsolineId[j] == i)
					{
						valid = true;
						break;
					}
					else if (strokeIsolineId[j] > i)
					{
						break;
					}
				}
				if (!valid)
				{
					m_isoline_finder->clearMarks(m_isoline, i);
				}
			}
			m_isoline = validIsolines;
			m_isoline_sequence.push_back(m_isoline);
		}
	}
}

int CMeshViewer::cutMesh()
{
	for (int id = 0; id < m_isoline_sequence.size(); id++)
	{
		m_isoline_finder->clearMarks(m_isoline_sequence[id]);
	}
	m_isoline_sequence.clear();

	m_isoline = m_isoline_finder->getIsoline(m_pMesh, 1.0 * m_selection / (m_isoline_set_size + 1));
	filterIsolines();
	for (int i = 0; i < m_isoline[0].size(); i++)
	{
		if (m_isoline[0].faces[i]->boundary_mark())
		{
			QMessageBox box(QMessageBox::Warning, "Intersecting boundaries", "Please redraw a stroke to avoid existing boundaries.", QMessageBox::Ok, NULL);
			box.exec();
			return 1;
		}
	}
	//m_isoline_finder->restoreMarks(m_pMesh, m_isoline);

	std::cout << "Selected isoline #" << m_selection << '\n' << std::endl;

	// mark boundary
	for (int i = 0; i < m_isoline[0].faces.size(); i++)
	{
		m_isoline[0].faces[i]->boundary_mark() = true;
	}
	for (int i = 0; i < m_isoline[0].edges.size(); i++)
	{
		m_isoline[0].edges[i]->bd_cross() = m_isoline[0].edges[i]->iso_cross();
	}
	m_boundary_faces.push_back(m_isoline[0].faces);

	m_showCut = true;
	// divide the mesh by coloring
	for (CVDMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++)
	{
		CViewerVertex *pV = *viter;
		pV->colormark() = false;
		if (m_updateMode)
		{
			pV->rgb() = pV->rgb_old();
			pV->part_label() = pV->part_label_prev();
		}
		else
		{
			pV->rgb_old() = pV->rgb();
			pV->part_label_prev() = pV->part_label();
		}
	}
	if (m_part_count == 1)
	{
		for (int j = 0; j < m_boundary_faces.back().size(); j++)
		{
			CViewerFace *pF = m_boundary_faces.back()[j];
			for (CVDMesh::FaceVertexIterator viter(pF); !viter.end(); viter++)
			{
				CViewerVertex *pV = *viter;
				pV->colormark() = true;
				if (pV->field() <= 1.0 * m_selection / (m_isoline_set_size + 1))
				{
					pV->rgb() = m_palette[0];
					pV->part_label() = 0;
				}
				else
				{
					pV->rgb() = m_palette[1];
					pV->part_label() = 1;
				}
			}
		}
		m_part_count++;
	}
	else
	{
		// ensure the new cut is far from other cuts
		if (m_strokes.curr_ends.start->part_label() != m_strokes.curr_ends.end->part_label())
		{
			QMessageBox box(QMessageBox::Warning, "Close boundaries", "Please redraw a stroke to avoid existing boundaries.", QMessageBox::Ok, NULL);
			box.exec();
			return 2;
		}
		int partToCut = m_strokes.curr_ends.start->part_label();
		// keep old colors for the regions not being cut
		for (CVDMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++)
		{
			CViewerVertex *pV = *viter;
			if (pV->part_label() != partToCut)
			{
				pV->colormark() = true;
			}
		}

		int highSide = 0;
		int lowSide = 0;
		std::vector<bool> sideLabel(m_part_count - 1, 0);
		for (int i = 0; i < m_part_count - 1; i++)
		{
			CViewerFace *pF = m_boundary_faces[i][0];
			CViewerHalfEdge *pH = m_pMesh->faceMostCcwHalfEdge(pF);
			CViewerVertex *pV = m_pMesh->halfedgeVertex(pH);
			if (pV->field() > 1.0 * m_selection / (m_isoline_set_size + 1))
			{
				highSide++;
				sideLabel[i] = 1;
			}
			else
			{
				lowSide++;
			}
		}
		if (highSide >= lowSide)
		{
			// deal with old boundaries
			for (int i = 0; i < m_part_count - 1; i++)
			{
				if (sideLabel[i])
				{
					for (int j = 0; j < m_boundary_faces[i].size(); j++)
					{
						CViewerFace *pF = m_boundary_faces[i][j];
						for (CVDMesh::FaceVertexIterator viter(pF); !viter.end(); viter++)
						{
							CViewerVertex *pV = *viter;
							pV->colormark() = true;
						}
					}
				}
				else
				{
					for (int j = 0; j < m_boundary_faces[i].size(); j++)
					{
						CViewerFace *pF = m_boundary_faces[i][j];
						for (CVDMesh::FaceVertexIterator viter(pF); !viter.end(); viter++)
						{
							CViewerVertex *pV = *viter;
							pV->colormark() = true;
							if (pV->part_label() == partToCut)
							{
								pV->rgb() = m_palette[m_part_count];
								pV->part_label() = m_part_count;
							}
						}
					}
				}
			}
			// deal with newly added boundary
			for (int j = 0; j < m_boundary_faces.back().size(); j++)
			{
				CViewerFace *pF = m_boundary_faces.back()[j];
				for (CVDMesh::FaceVertexIterator viter(pF); !viter.end(); viter++)
				{
					CViewerVertex *pV = *viter;
					pV->colormark() = true;
					if (pV->field() <= 1.0 * m_selection / (m_isoline_set_size + 1))
					{
						pV->rgb() = m_palette[m_part_count];
						pV->part_label() = m_part_count;
					}
				}
			}
		}
		else
		{
			// deal with old boundaries
			for (int i = 0; i < m_part_count - 1; i++)
			{
				if (!sideLabel[i])
				{
					for (int j = 0; j < m_boundary_faces[i].size(); j++)
					{
						CViewerFace *pF = m_boundary_faces[i][j];
						for (CVDMesh::FaceVertexIterator viter(pF); !viter.end(); viter++)
						{
							CViewerVertex *pV = *viter;
							pV->colormark() = true;
						}
					}
				}
				else
				{
					for (int j = 0; j < m_boundary_faces[i].size(); j++)
					{
						CViewerFace *pF = m_boundary_faces[i][j];
						for (CVDMesh::FaceVertexIterator viter(pF); !viter.end(); viter++)
						{
							CViewerVertex *pV = *viter;
							pV->colormark() = true;
							if (pV->part_label() == partToCut)
							{
								pV->rgb() = m_palette[m_part_count];
								pV->part_label() = m_part_count;
							}
						}
					}
				}
			}
			// deal with newly added boundary
			for (int j = 0; j < m_boundary_faces.back().size(); j++)
			{
				CViewerFace *pF = m_boundary_faces.back()[j];
				for (CVDMesh::FaceVertexIterator viter(pF); !viter.end(); viter++)
				{
					CViewerVertex *pV = *viter;
					pV->colormark() = true;
					if (pV->field() > 1.0 * m_selection / (m_isoline_set_size + 1))
					{
						pV->rgb() = m_palette[m_part_count];
						pV->part_label() = m_part_count;
					}
				}
			}
		}
		m_part_count++;
	}

	bool ending = true;
	do
	{
		ending = true;
		for (CVDMesh::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); eiter++)
		{
			CViewerEdge *pE = *eiter;
			CViewerVertex *pV = m_pMesh->edgeVertex1(pE);
			CViewerVertex *pW = m_pMesh->edgeVertex2(pE);
			if (pV->colormark() && !pW->colormark())
			{
				ending = false;
				dyeNeighborhood(pV);
				break;
			}
			else if (!pV->colormark() && pW->colormark())
			{
				ending = false;
				dyeNeighborhood(pW);
				break;
			}
			else
			{
				continue;
			}
		}
	} while (!ending);

	return 0;
}

void CMeshViewer::dyeNeighborhood(CViewerVertex *pV)
{
	for (CVDMesh::VertexVertexIterator witer(pV); !witer.end(); witer++)
	{
		CViewerVertex *pW = *witer;
		if (pW->colormark())
		{
			continue;
		}
		else
		{
			pW->colormark() = true;
			pW->rgb() = pV->rgb();
			pW->part_label() = pV->part_label();
			dyeNeighborhood(pW);
		}
	}
}

FL CMeshViewer::findPath(CViewerFace *start, CViewerFace *end)
{
	for (CVDMesh::MeshFaceIterator fiter(m_pMesh); !fiter.end(); fiter++)
	{
		CViewerFace *pF = *fiter;
		pF->layer() = 0;
		pF->parent() = nullptr;
	}

	std::queue<CViewerFace*> faceQ;
	start->layer() = 1;
	faceQ.push(start);
	bool ending = false;
	while (!ending && !faceQ.empty())
	{
		CViewerFace *pF = faceQ.front();
		for (CVDMesh::FaceVertexIterator viter(pF); !viter.end(); viter++)
		{
			CViewerVertex *pV = *viter;
			for (CVDMesh::VertexFaceIterator giter(pV); !giter.end(); giter++)
			{
				CViewerFace *pG = *giter;
				if (pG->layer() == 0)
				{
					pG->layer() = pF->layer() + 1;
					pG->parent() = pF;
					faceQ.push(pG);
				}
				if (pG == end)
				{
					ending = true;
					break;
				}
			}
			if (ending)
			{
				break;
			}
		}
		faceQ.pop();
	}

	FL faceV;
	faceV.push_back(end);
	CViewerFace *f = end->parent();
	while (f != start)
	{
		f->picked() = true;
		faceV.push_back(f);
		f = f->parent();
	}
	faceV.push_back(start);
	std::reverse(faceV.begin(), faceV.end());
	return faceV;
}

void CMeshViewer::regionGrow()
{
	if (!m_strokes.all_end_faces.empty())
	{
		FL marked, candidates;

		// initialize
		for (CVDMesh::MeshFaceIterator fiter(m_pMesh); !fiter.end(); fiter++)
		{
			CViewerFace *pF = *fiter;
			pF->patch_label() = -1;
			pF->isCandidate() = false;
			pF->marked_neighbor() = nullptr;
			pF->metric() = INFI;
		}
		for (int i = 0; i < m_strokes.all_end_faces.size(); i++)
		{
			CViewerFace *f0 = m_strokes.all_end_faces[i].start;
			CViewerFace *f1 = m_strokes.all_end_faces[i].end;
			f0->patch_label() = 0;
			f0->rgb() = m_palette[0];
			f1->patch_label() = 1;
			f1->rgb() = m_palette[1];
			marked.push_back(f0);
			marked.push_back(f1);
		}
		for (int i = 0; i < marked.size(); i++)
		{
			CViewerFace *pF = marked[i];
			for (CVDMesh::FaceHalfedgeIterator hiter(pF); !hiter.end(); hiter++)
			{
				CViewerEdge *pE = m_pMesh->halfedgeEdge(*hiter);
				CViewerFace *pG = m_pMesh->eOtherFace(pE, pF);
				// only deal with unmarked faces
				if (pG->patch_label() == -1)
				{
					double cost_ij = cost(pF, pG);
					if (cost_ij < pG->metric())
					{
						pG->metric() = cost_ij;
						pG->marked_neighbor() = pF;
						if (!pG->isCandidate())
						{
							pG->isCandidate() = true;
							candidates.push_back(pG);
						}
					}
				}
			}
		}

		// iteratively absorb an unmarked neighboring face with min cost
		while (!candidates.empty())
		{
			// find the face with min cost and absorb it
			double minCost = INFI;
			CViewerFace *pF_min = nullptr;
			int minID = 0;
			for (int i = 0; i < candidates.size(); i++)
			{
				CViewerFace *pF = candidates[i];
				if (pF->metric() < minCost)
				{
					minCost = pF->metric();
					pF_min = pF;
					minID = i;
				}
			}
			pF_min->patch_label() = pF_min->marked_neighbor()->patch_label();
			pF_min->rgb() = pF_min->marked_neighbor()->rgb();

			// update candidate list
			for (CVDMesh::FaceHalfedgeIterator hiter(pF_min); !hiter.end(); hiter++)
			{
				CViewerEdge *pE = m_pMesh->halfedgeEdge(*hiter);
				CViewerFace *pG = m_pMesh->eOtherFace(pE, pF_min);
				// only deal with unmarked faces
				if (pG->patch_label() == -1)
				{
					double cost_ij = cost(pF_min, pG);
					if (cost_ij < pG->metric())
					{
						pG->metric() = cost_ij;
						pG->marked_neighbor() = pF_min;
						if (!pG->isCandidate())
						{
							pG->isCandidate() = true;
							candidates.push_back(pG);
						}
					}
				}
			}
			pF_min->isCandidate() = false;
			candidates.erase(candidates.begin() + minID);
		}
	}
}

double CMeshViewer::cost(CViewerFace *fi, CViewerFace *fj)
{
	// here fi is the marked face and fj unmarked
	CPoint ni, nj;
	if (m_doNormalFilter)
	{
		ni = fi->normal_filtered();
		nj = fj->normal_filtered();
	}
	else
	{
		ni = fi->normal();
		nj = fj->normal();
	}
	CPoint n_ij = ni - nj;
	CPoint e_ij = m_pMesh->faceCenter(fi) - m_pMesh->faceCenter(fj);
	double theta_ij = (ni * e_ij) / e_ij.norm();
	double cost_ij = (1 + fabs(theta_ij)) * n_ij.norm();
	return cost_ij;
}

void CMeshViewer::normalFilter()
{
	for (CVDMesh::MeshFaceIterator fiter(m_pMesh); !fiter.end(); fiter++)
	{
		CViewerFace *fi = *fiter;
		// weighted sum of face normals in 2-ring neighborhood of fi
		for (CVDMesh::FaceVertexIterator viter(fi); !viter.end(); viter++)
		{
			CViewerVertex *pV = *viter;
			for (CVDMesh::VertexFaceIterator giter(pV); !giter.end(); giter++)
			{
				CViewerFace *fj = *giter;
				if (fj != fi)
				{
					fi->normal_filtered() += fj->normal() * GaussianWeight(fi, fj);
				}
			}
		}
		fi->normal_filtered() /= fi->normal_filtered().norm();
	}
}

double CMeshViewer::GaussianWeight(CViewerFace *fi, CViewerFace *fj)
{
	CPoint n_ij = fi->normal() - fj->normal();
	double sigma = 0.35;
	return exp(-n_ij.norm() * n_ij.norm() / (2 * sigma * sigma));
}