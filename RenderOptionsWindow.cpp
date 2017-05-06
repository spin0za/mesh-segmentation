#include <QtWidgets>

#include "RenderOptionsWindow.h"

RenderOptionsWindow::RenderOptionsWindow(CMeshViewer * pViewer, QWidget *parent)
: QWidget(parent)
{
	m_pViewer = pViewer;

	QGridLayout *grid = new QGridLayout;

	grid->addWidget(createTextureGroup(), 0, 0);
	grid->addWidget(createNormalGroup(), 1, 0);
	grid->addWidget(createGeometryGroup(), 2, 0);
	grid->addWidget(createPickGroup(), 3, 0);


	grid->addWidget(createShowOptions(), 0, 1);
	grid->addWidget(createManipulationGroup(), 1, 1);
	grid->addWidget(createPolygonGroup(), 2, 1);
	grid->addWidget(createBrushGroup(), 3, 1);
	setLayout(grid);

	setWindowTitle(tr("Render Options"));
	resize(480, 320);
}

RenderOptionsWindow::RenderOptionsWindow(QWidget *parent)
: QWidget(parent)
{
	QGridLayout *grid = new QGridLayout;

	grid->addWidget(createTextureGroup(), 0, 0);
	grid->addWidget(createNormalGroup(), 1, 0);
	grid->addWidget(createGeometryGroup(), 2, 0);
	

	grid->addWidget(createShowOptions(), 0, 1);
	grid->addWidget(createManipulationGroup(), 1, 1);
	grid->addWidget(createPolygonGroup(), 2, 1);
	setLayout(grid);

	setWindowTitle(tr("Render Options"));
	resize(480, 320);
}

QGroupBox *RenderOptionsWindow::createTextureGroup()
{
	texture_group = new QGroupBox(tr("Texture Map Mode"));
	texture_group->setCheckable(true);

	texture_replace    = new QRadioButton(tr("&Replace"));
	texture_modulate   = new QRadioButton(tr("&Modulation"));
	texture_none	   = new QRadioButton(tr("&None"));

	if (m_pViewer != NULL && m_pViewer->get_texture_image() != NULL)
	{
		texture_group->setCheckable(true);

		switch ( m_pViewer->get_texture_mode())
		{
		case TEXTURE_MODE::MODULATE:
			texture_modulate->setChecked(true);
			break;
		case TEXTURE_MODE::REPLACE:
			texture_replace->setChecked(true);
			break;
		case TEXTURE_MODE::NONE:
			texture_none->setChecked(true);
			break;
		}
	}
	else
	{
		texture_group->setChecked(false);
	}

	QVBoxLayout *vbox = new QVBoxLayout;
	vbox->addWidget(texture_replace);
	vbox->addWidget(texture_modulate);
	vbox->addWidget(texture_none);
	vbox->addStretch(1);
	texture_group->setLayout(vbox);

	connect(texture_group, SIGNAL(clicked(bool)), this, SLOT(textureGroupClicked(bool)));
	connect(texture_replace, SIGNAL(clicked(bool)), this, SLOT(textureclicked()));
	connect(texture_modulate, SIGNAL(clicked(bool)), this, SLOT(textureclicked()));
	connect(texture_none, SIGNAL(clicked(bool)), this, SLOT(textureclicked()));

	return texture_group;
}

void RenderOptionsWindow::textureGroupClicked(bool b)
{
	if (b )
	{
		if (m_pViewer != NULL && m_pViewer->get_texture_image() != NULL)
		{

			texture_group->setChecked(true);
			switch (m_pViewer->get_texture_mode())
			{
			case TEXTURE_MODE::MODULATE:
				texture_modulate->setChecked(true);
				break;
			case TEXTURE_MODE::REPLACE:
				texture_replace->setChecked(true);
				break;
			case TEXTURE_MODE::NONE:
				texture_none->setChecked(true);
				break;
			}
		}
		else
		{
			texture_group->setChecked(false);
		}
	}
	else
	{
		texture_group->setChecked(false);
		m_pViewer->set_texture_mode(TEXTURE_MODE::NONE);
	}
}

void RenderOptionsWindow::textureclicked()
{
	if (texture_replace->isChecked())
	{
		std::cout << "replace" << std::endl;
		m_pViewer->set_texture_mode(TEXTURE_MODE::REPLACE);
	}
	if (texture_modulate->isChecked())
	{
		std::cout << "modulate" << std::endl;
		m_pViewer->set_texture_mode(TEXTURE_MODE::MODULATE);
	}
	if (texture_none->isChecked())
	{
		std::cout << "texture disabled" << std::endl;
		m_pViewer->set_texture_mode(TEXTURE_MODE::NONE);
	}
}

QGroupBox *RenderOptionsWindow::createNormalGroup()
{
	normal_group = new QGroupBox(tr("&Normal Mode"));
	

	normal_vertex = new QRadioButton(tr("&Vertex normal"));
	normal_face   = new QRadioButton(tr("&Face normal"));
	if (!m_pViewer)
	{
		normal_vertex->setChecked(true);
	}
	else
	{
		switch (m_pViewer->get_normal_mode())
		{
		case NORMAL_MODE::FACE:
			normal_face->setChecked(true);
			break;
		case NORMAL_MODE::VERTEX:
			normal_vertex->setChecked(true);
			break;
		}
	}
	QVBoxLayout *vbox = new QVBoxLayout;
	vbox->addWidget(normal_vertex);
	vbox->addWidget(normal_face);
	vbox->addStretch(1);
	normal_group->setLayout(vbox);

	connect(normal_vertex, SIGNAL(clicked(bool)), this, SLOT(normalclicked()));
	connect(normal_face,   SIGNAL(clicked(bool)), this, SLOT(normalclicked()));

	return normal_group;
}

void RenderOptionsWindow::normalclicked()
{
	if (normal_vertex->isChecked())
	{
		std::cout << "vertex" << std::endl;
		m_pViewer->set_normal_mode( NORMAL_MODE::VERTEX );

	}
	if (normal_face->isChecked())
	{
		std::cout << "face" << std::endl;
		m_pViewer->set_normal_mode(NORMAL_MODE::FACE);
	}

}


QGroupBox *RenderOptionsWindow::createGeometryGroup()
{
	geometry_group = new QGroupBox(tr("&Geometry Mode"));
	
	geometry_space = new QRadioButton(tr("&Spacial Domain"));
	geometry_plane = new QRadioButton(tr("&Parameter Domain"));
	
	if (!m_pViewer)
	{
		geometry_space->setChecked(true);
	}
	else
	{
		switch (m_pViewer->get_geometry_mode())
		{
		case GEOMETRY_MODE::GEOMETRY:
			geometry_space->setChecked(true);
			break;
		case GEOMETRY_MODE::PARAMETER:
			geometry_plane->setChecked(true);
			break;
		}
	}


	QVBoxLayout *vbox = new QVBoxLayout;
	vbox->addWidget(geometry_space);
	vbox->addWidget(geometry_plane);
	vbox->addStretch(1);
	geometry_group->setLayout(vbox);

	connect(geometry_space, SIGNAL(clicked(bool)), this, SLOT(geometryclicked()));
	connect(geometry_plane, SIGNAL(clicked(bool)), this, SLOT(geometryclicked()));

	return geometry_group;
}

void RenderOptionsWindow::geometryclicked()
{
	if (geometry_space->isChecked())
	{
		std::cout << "Geometry Space" << std::endl;
		m_pViewer->set_geometry_mode(GEOMETRY_MODE::GEOMETRY); 
	}
	if (geometry_plane->isChecked())
	{
		std::cout << "Geometry Plane" << std::endl;
		m_pViewer->set_geometry_mode(GEOMETRY_MODE::PARAMETER);
	}
}

QGroupBox *RenderOptionsWindow::createPolygonGroup()
{
	polygon_group = new QGroupBox(tr("Polygon Mode"));

	polygon_line = new QRadioButton(tr("&Line"));
	polygon_fill = new QRadioButton(tr("&Fill"));

	if (!m_pViewer)
	{
		polygon_fill->setChecked(true);
	}
	else
	{
		switch (m_pViewer->get_polygon_mode())
		{
		case POLYGON_MODE::FILL:
			polygon_fill->setChecked(true);
			break;
		case POLYGON_MODE::LINE:
			polygon_line->setChecked(true);
			break;
		}
	}


	QVBoxLayout *vbox = new QVBoxLayout;
	vbox->addWidget(polygon_fill);
	vbox->addWidget(polygon_line);
	vbox->addStretch(1);
	polygon_group->setLayout(vbox);

	connect(polygon_fill, SIGNAL(clicked(bool)), this, SLOT(polygonclicked()));
	connect(polygon_line, SIGNAL(clicked(bool)), this, SLOT(polygonclicked()));

	return polygon_group;
}

void RenderOptionsWindow::polygonclicked()
{

	if (polygon_fill->isChecked())
	{
		std::cout << "polygon fill" << std::endl;
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		m_pViewer->set_polygon_mode(POLYGON_MODE::FILL);
	}
	if (polygon_line->isChecked())
	{
		std::cout << "polygon line" << std::endl;
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		m_pViewer->set_polygon_mode(POLYGON_MODE::LINE);
	}
}

QGroupBox *RenderOptionsWindow::createManipulationGroup()
{
	manipulation_group = new QGroupBox(tr("Manipulation Mode"));

	manipulation_eye = new QRadioButton(tr("&Camera"));
	manipulation_obj = new QRadioButton(tr("&Object"));

	if (!m_pViewer)
	{
		manipulation_obj->setChecked(true);
	}
	else
	{
		switch (m_pViewer->get_manipulation_mode())
		{
		case MANIPULATION_MODE::EYE:
			manipulation_eye->setChecked(true);
			break;
		case MANIPULATION_MODE::OBJECT:
			manipulation_obj->setChecked(true);
			break;
		}
	}

	QVBoxLayout *vbox = new QVBoxLayout;
	vbox->addWidget(manipulation_eye);
	vbox->addWidget(manipulation_obj);
	vbox->addStretch(1);
	manipulation_group->setLayout(vbox);

	connect(manipulation_eye, SIGNAL(clicked(bool)), this, SLOT(manipulationclicked()));
	connect(manipulation_obj, SIGNAL(clicked(bool)), this, SLOT(manipulationclicked()));
	return manipulation_group;
}


QGroupBox *RenderOptionsWindow::createPickGroup()
{
	pick_group = new QGroupBox(tr("Picking Mode"));

	pick_vertex = new QRadioButton(tr("&Vertex"));
	pick_face   = new QRadioButton(tr("&Face"));

	if (!m_pViewer)
	{
		pick_vertex->setChecked(true);
	}
	else
	{
		switch (m_pViewer->get_pick_mode())
		{
		case PICK_MODE::VERTEX:
			pick_vertex->setChecked(true);
			break;
		case PICK_MODE::FACE:
			pick_face->setChecked(true);
			break;
		}
	}

	QVBoxLayout *vbox = new QVBoxLayout;
	vbox->addWidget(pick_vertex);
	vbox->addWidget(pick_face);
	vbox->addStretch(1);
	pick_group->setLayout(vbox);

	connect(pick_vertex, SIGNAL(clicked(bool)), this, SLOT(pickclicked()));
	connect(pick_face, SIGNAL(clicked(bool)), this, SLOT(pickclicked()));
	return pick_group;
}

QGroupBox *RenderOptionsWindow::createBrushGroup()
{
	brush_group = new QGroupBox(tr("Brush Mode"));

	brush_part = new QRadioButton(tr("&Part"));
	brush_patch = new QRadioButton(tr("&Patch"));

	if (!m_pViewer)
	{
		brush_part->setChecked(true);
	}
	else
	{
		switch (m_pViewer->get_brush_mode())
		{
		case BRUSH_MODE::PART:
			brush_part->setChecked(true);
			break;
		case BRUSH_MODE::PATCH:
			brush_patch->setChecked(true);
			break;
		}
	}

	QVBoxLayout *vbox = new QVBoxLayout;
	vbox->addWidget(brush_part);
	vbox->addWidget(brush_patch);
	vbox->addStretch(1);
	brush_group->setLayout(vbox);

	connect(brush_part, SIGNAL(clicked(bool)), this, SLOT(brushclicked()));
	connect(brush_patch, SIGNAL(clicked(bool)), this, SLOT(brushclicked()));
	return brush_group;
}

void RenderOptionsWindow::manipulationclicked()
{
	if (manipulation_eye->isChecked())
	{
		m_pViewer->set_manipulation_mode(MANIPULATION_MODE::EYE);
		std::cout << "Manipulation Eye" << std::endl;

	}
	if (manipulation_obj->isChecked())
	{
		m_pViewer->set_manipulation_mode(MANIPULATION_MODE::OBJECT);
		std::cout << "Manipulation Object" << std::endl;
	}
}

void RenderOptionsWindow::pickclicked()
{
	if (pick_vertex->isChecked())
	{
		m_pViewer->set_pick_mode(PICK_MODE::VERTEX);
		std::cout << "Pick Vertex" << std::endl;
	}
	if (pick_face->isChecked())
	{
		m_pViewer->set_pick_mode(PICK_MODE::FACE);
		std::cout << "Pick Face" << std::endl;
	}
}

void RenderOptionsWindow::brushclicked()
{
	if (brush_part->isChecked())
	{
		m_pViewer->set_brush_mode(BRUSH_MODE::PART);
		m_pViewer->set_face_rgb_flag(false);
		std::cout << "Part Brush" << std::endl;
	}
	if (brush_patch->isChecked())
	{
		m_pViewer->set_brush_mode(BRUSH_MODE::PATCH);
		m_pViewer->set_face_rgb_flag(true);
		m_pViewer->set_show_cut_flag(false);
		std::cout << "Patch Brush" << std::endl;
	}
}

QGroupBox *RenderOptionsWindow::createShowOptions()
{
	QGroupBox *groupBox = new QGroupBox(tr("Show Options"));
	//groupBox->setChecked(false);


	show_edge = new QCheckBox(tr("Show Edges"));
	show_boundary   = new QCheckBox(tr("Show Boundary"));
	show_backface   = new QCheckBox(tr("Show Backface"));
	show_sharp_edge = new QCheckBox(tr("Show Sharp Edges"));
	normal_filter	= new QCheckBox(tr("Normal Filter"));
	normal_reverse  = new QCheckBox(tr("Normal Reverse"));
	show_vertex_rgb = new QCheckBox(tr("Show Vertex Rgb"));
	show_picked = new QCheckBox(tr("Show Picked Vertex/Face"));
	slow_motion = new QCheckBox(tr("Show Intermediate Results"));


	if (!m_pViewer)
	{
		show_edge->setChecked(false);
		show_boundary->setChecked(false);
		show_backface->setChecked(false);
		normal_filter->setChecked(false);
		normal_reverse->setChecked(false);
		show_sharp_edge->setChecked(false);
		show_vertex_rgb->setChecked(false);
		show_picked->setChecked(false);
		slow_motion->setChecked(false);
	}
	else
	{
		show_edge->setChecked( m_pViewer->get_edge_flag() );
		show_boundary->setChecked(m_pViewer->get_boundary_flag());
		show_backface->setChecked(m_pViewer->get_backface_flag());
		normal_filter->setChecked(m_pViewer->get_filter_normal_flag());
		normal_reverse->setChecked(m_pViewer->get_reverse_normal_flag());
		show_sharp_edge->setChecked(m_pViewer->get_sharp_edge_flag());
		show_vertex_rgb->setChecked(m_pViewer->get_vertex_rgb_flag());
		show_picked->setChecked(m_pViewer->get_picked_flag());
		slow_motion->setChecked(m_pViewer->get_slow_motion_flag());
	}

	QVBoxLayout *vbox = new QVBoxLayout;
	vbox->addWidget(show_edge);
	vbox->addWidget(show_sharp_edge); 
	vbox->addWidget(show_boundary);
	vbox->addWidget(show_backface);
	vbox->addWidget(normal_filter);
	vbox->addWidget(normal_reverse);
	vbox->addWidget(show_vertex_rgb);
	vbox->addWidget(show_picked);
	vbox->addWidget(slow_motion);

	vbox->addStretch(1);
	groupBox->setLayout(vbox);


	connect(show_edge,	   SIGNAL(clicked(bool)), this, SLOT(ShowEdgeClicked(bool)));
	connect(show_sharp_edge, SIGNAL(clicked(bool)), this, SLOT(ShowSharpEdgeClicked(bool)));
	connect(show_boundary, SIGNAL(clicked(bool)), this, SLOT(ShowBoundaryClicked(bool)));
	connect(show_backface, SIGNAL(clicked(bool)), this, SLOT(ShowBackFaceClicked(bool)));
	connect(normal_filter, SIGNAL(clicked(bool)), this, SLOT(NormalFilterClicked(bool)));
	connect(normal_reverse,SIGNAL(clicked(bool)), this, SLOT(NormalReverseClicked(bool)));
	connect(show_vertex_rgb, SIGNAL(clicked(bool)), this, SLOT(ShowVertexRgbClicked(bool)));
	connect(show_picked, SIGNAL(clicked(bool)), this, SLOT(ShowPickedClicked(bool)));
	connect(slow_motion, SIGNAL(clicked(bool)), this, SLOT(SlowMotionClicked(bool)));

	return groupBox;
}

void RenderOptionsWindow::NormalFilterClicked(bool b)
{
	if (b)
	{
		std::cout << "Enabled Normal Filtering" << std::endl;
		m_pViewer->set_filter_normal_flag(true);
	}
	else
	{
		std::cout << "Disabled Normal Filtering" << std::endl;
		m_pViewer->set_filter_normal_flag(false);
	}
}

void RenderOptionsWindow::NormalReverseClicked(bool b)
{
	if (b)
	{
		std::cout << "Normal Reverse" << std::endl;
		m_pViewer->set_reverse_normal_flag(true);
	}
	else
	{
		std::cout << "Original Normal" << std::endl;
		m_pViewer->set_reverse_normal_flag(false);
	}
}

void RenderOptionsWindow::ShowBackFaceClicked(bool b)
{
	if (b)
	{
		std::cout << "Show BackFace" << std::endl;
		m_pViewer->set_backface_flag(true);
	}
	else
	{
		std::cout << "Hide BaceFace" << std::endl;
		m_pViewer->set_backface_flag(false);
	}
}

void RenderOptionsWindow::ShowPickedClicked(bool b)
{
	if (b)
	{
		std::cout << "Show Picked Vertex/Face" << std::endl;
		m_pViewer->set_picked_flag(true);
	}
	else
	{
		std::cout << "Hide Picked Vertex/Face" << std::endl;
		m_pViewer->set_picked_flag(false);
	}
}


void RenderOptionsWindow::SlowMotionClicked(bool b)
{
	if (b)
	{
		std::cout << "Show intermediate results" << std::endl;
		m_pViewer->set_slow_motion_flag(true);
	}
	else
	{
		std::cout << "Skip intermediate results" << std::endl;
		m_pViewer->set_slow_motion_flag(false);
	}
}

void RenderOptionsWindow::ShowBoundaryClicked(bool b)
{
	if (b)
	{
		std::cout << "Show Boundary" << std::endl;
		m_pViewer->set_boundary_flag(true);
	}
	else
	{
		std::cout << "Hide Boundary" << std::endl;
		m_pViewer->set_boundary_flag(false);
	}
}

void RenderOptionsWindow::ShowEdgeClicked(bool b)
{
	if (b)
	{
		std::cout << "Show Edge" << std::endl;
		m_pViewer->set_edge_flag(true);
	}
	else
	{
		std::cout << "Hide Edge" << std::endl;
		m_pViewer->set_edge_flag(false);
	}
}

void RenderOptionsWindow::ShowSharpEdgeClicked(bool b)
{
	if (b)
	{
		std::cout << "Show Sharp Edge" << std::endl;
		m_pViewer->set_sharp_edge_flag(true);
	}
	else
	{
		std::cout << "Hide Sharp Edge" << std::endl;
		m_pViewer->set_sharp_edge_flag(false);
	}
}

void RenderOptionsWindow::ShowVertexRgbClicked(bool b)
{
	if (b)
	{
		std::cout << "Show Vertex Rgb" << std::endl;
		m_pViewer->set_vertex_rgb_flag(true);
	}
	else
	{
		std::cout << "Hide Vertex Rgb" << std::endl;
		m_pViewer->set_vertex_rgb_flag(false);
	}
}

QGroupBox *RenderOptionsWindow::createNonExclusiveGroup()
{
	QGroupBox *groupBox = new QGroupBox(tr("Non-Exclusive Checkboxes"));
	groupBox->setFlat(true);

	QCheckBox *checkBox1 = new QCheckBox(tr("&Checkbox 1"));
	QCheckBox *checkBox2 = new QCheckBox(tr("C&heckbox 2"));
	checkBox2->setChecked(true);
	QCheckBox *tristateBox = new QCheckBox(tr("Tri-&state button"));
	tristateBox->setTristate(true);
	tristateBox->setCheckState(Qt::PartiallyChecked);

	QVBoxLayout *vbox = new QVBoxLayout;
	vbox->addWidget(checkBox1);
	vbox->addWidget(checkBox2);
	vbox->addWidget(tristateBox);
	vbox->addStretch(1);
	groupBox->setLayout(vbox);

	return groupBox;
}



QGroupBox *RenderOptionsWindow::createPushButtonGroup()
{
	QGroupBox *groupBox = new QGroupBox(tr("&Push Buttons"));
	groupBox->setCheckable(true);
	groupBox->setChecked(true);

	QPushButton *pushButton = new QPushButton(tr("&Normal Button"));
	QPushButton *toggleButton = new QPushButton(tr("&Toggle Button"));
	toggleButton->setCheckable(true);
	toggleButton->setChecked(true);
	QPushButton *flatButton = new QPushButton(tr("&Flat Button"));
	flatButton->setFlat(true);

	QPushButton *popupButton = new QPushButton(tr("Pop&up Button"));
	QMenu *menu = new QMenu(this);
	menu->addAction(tr("&First Item"));
	menu->addAction(tr("&Second Item"));
	menu->addAction(tr("&Third Item"));
	menu->addAction(tr("F&ourth Item"));
	popupButton->setMenu(menu);

	QAction *newAction = menu->addAction(tr("Submenu"));
	QMenu *subMenu = new QMenu(tr("Popup Submenu"));
	subMenu->addAction(tr("Item 1"));
	subMenu->addAction(tr("Item 2"));
	subMenu->addAction(tr("Item 3"));
	newAction->setMenu(subMenu);

	QVBoxLayout *vbox = new QVBoxLayout;
	vbox->addWidget(pushButton);
	vbox->addWidget(toggleButton);
	vbox->addWidget(flatButton);
	vbox->addWidget(popupButton);
	vbox->addStretch(1);
	groupBox->setLayout(vbox);

	return groupBox;
}