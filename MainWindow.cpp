#include "MainWindow.h"

MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	m_pViewer = new CMeshViewer();
	setCentralWidget(m_pViewer);
	createActions();
	createToolbar();

	setAcceptDrops(true);
		
	updateCursor = new QCursor(QPixmap(":/icons/images/cross48.png"), -1, -1);
}

MainWindow::~MainWindow()
{
}

void MainWindow::createActions()
{
	openMeshAction = new QAction(tr("&Open Mesh"), this);
	openMeshAction->setIcon(QIcon(":/icons/images/open.png"));
	openMeshAction->setShortcut(QKeySequence::Open);
	openMeshAction->setStatusTip(tr("Open a mesh file"));
	connect(openMeshAction, SIGNAL(triggered()), m_pViewer, SLOT(openMesh()));

	openTextureAction = new QAction(tr("&Open Texture"), this);
	openTextureAction->setIcon(QIcon(":/icons/images/texture.png"));
	//openTextureAction->setShortcut(QKeySequence::Open);
	openTextureAction->setStatusTip(tr("Open a texture file"));
	connect(openTextureAction, SIGNAL(triggered()), m_pViewer, SLOT(openTexture()));


	viewPoints = new QAction(tr("&Points"), this);
	viewPoints->setIcon(QIcon(":/icons/images/points.png"));
	viewPoints->setText(tr("Draw Points"));
	viewPoints->setStatusTip(tr("Points"));
	viewPoints->setCheckable(true);
	viewPoints->setChecked(false);
	connect(viewPoints, SIGNAL(triggered()), m_pViewer, SLOT(showPoints()));

	renderOptions = new QAction(tr("&Render"), this);
	renderOptions->setIcon(QIcon(":/icons/images/render.png"));
	renderOptions->setText(tr("Set render options"));
	renderOptions->setStatusTip(tr("Render"));
	connect(renderOptions, SIGNAL(triggered()), this, SLOT(setRenderOptions()));

	viewWireframe = new QAction(tr("&Wireframe"), this);
	viewWireframe->setIcon(QIcon(":/icons/images/wireframe.png"));
	viewWireframe->setText(tr("Draw Wireframe"));
	viewWireframe->setStatusTip(tr("Wireframe"));
	viewWireframe->setCheckable(true);
	viewWireframe->setChecked(false);
	connect(viewWireframe, SIGNAL(triggered()), m_pViewer, SLOT(showWireframe()));

	viewFlatlines = new QAction(tr("Flat&lines"), this);
	viewFlatlines->setIcon(QIcon(":/icons/images/flatlines.png"));
	viewFlatlines->setText(tr("Draw Flatlines"));
	viewFlatlines->setStatusTip(tr("Flatlines"));
	viewFlatlines->setCheckable(true);
	viewFlatlines->setChecked(false);
	connect(viewFlatlines, SIGNAL(triggered()), m_pViewer, SLOT(showFlatlines()));

	viewFlat = new QAction(tr("&Flat"), this);
	viewFlat->setIcon(QIcon(":/icons/images/flat.png"));
	viewFlat->setText(tr("Draw Flat"));
	viewFlat->setStatusTip(tr("Flat"));
	viewFlat->setCheckable(true);
	viewFlat->setChecked(false);
	connect(viewFlat, SIGNAL(triggered()), m_pViewer, SLOT(showFlat()));

	viewSmooth = new QAction(tr("&Smooth"), this);
	viewSmooth->setIcon(QIcon(":/icons/images/smooth.png"));
	viewSmooth->setText(tr("Draw Smooth"));
	viewSmooth->setStatusTip(tr("Smooth"));
	viewSmooth->setCheckable(true);
	viewSmooth->setChecked(false);
	connect(viewSmooth, SIGNAL(triggered()), m_pViewer, SLOT(showSmooth()));

	viewBoundary = new checkableAction(this);
	viewBoundary->setIcon(QIcon(":/icons/images/boundary.png"));
	viewBoundary->setText(tr("Draw Boundary"));
	viewBoundary->setStatusTip(tr("Boundary"));
	viewBoundary->setCheckable(true);
	viewBoundary->setChecked(false);
	connect(viewBoundary, SIGNAL(actionCheck()), m_pViewer, SLOT(showBoundary()));
	connect(viewBoundary, SIGNAL(actionUncheck()), m_pViewer, SLOT(showBoundary()));

	viewTexture = new QAction(tr("&Texture"), this);
	viewTexture->setIcon(QIcon(":/icons/images/bunny-alpha0.png"));
	viewTexture->setText(tr("Draw Texture"));
	viewTexture->setStatusTip(tr("Texture"));
	viewTexture->setCheckable(true);
	viewTexture->setChecked(false);
	connect(viewTexture, SIGNAL(triggered()), m_pViewer, SLOT(showTexture()));

	takeSnapshot = new QAction(tr("&Snapshot"), this);
	takeSnapshot->setIcon(QIcon(":/icons/images/camera.png"));
	takeSnapshot->setText(tr("Take snapshot"));
	takeSnapshot->setStatusTip(tr("Snapshot"));
	connect(takeSnapshot, SIGNAL(triggered()), m_pViewer, SLOT(saveImage()));

	selectAction = new checkableAction(this);
	selectAction->setIcon(QIcon(":/icons/images/select.png"));
	selectAction->setText(tr("Select Vertices"));
	selectAction->setStatusTip(tr("Select Vertices"));
	selectAction->setCheckable(true);
	selectAction->setChecked(false);
	connect(selectAction, SIGNAL(actionCheck()), m_pViewer, SLOT(enterSelectionMode()));
	connect(selectAction, SIGNAL(actionUncheck()), m_pViewer, SLOT(quitSelectionMode()));

	strokeBrush = new checkableAction(this);
	strokeBrush->setIcon(QIcon(":/icons/images/brush.png"));
	strokeBrush->setText(tr("Draw a stroke"));
	strokeBrush->setStatusTip(tr("Draw a stroke"));
	strokeBrush->setCheckable(true);
	strokeBrush->setChecked(false);
	connect(strokeBrush, SIGNAL(actionCheck()), m_pViewer, SLOT(enterBrushMode()));
	connect(strokeBrush, SIGNAL(actionCheck()), this, SLOT(enterBrush()));
	connect(strokeBrush, SIGNAL(actionUncheck()), m_pViewer, SLOT(quitBrushMode()));
	connect(strokeBrush, SIGNAL(actionUncheck()), this, SLOT(quitBrush()));

	connect(strokeBrush, SIGNAL(actionCheck()), this, SLOT(quitSelection()));
	connect(strokeBrush, SIGNAL(actionCheck()), selectAction, SIGNAL(actionUncheck()));
	connect(selectAction, SIGNAL(actionCheck()), this, SLOT(quitBrush()));
	connect(selectAction, SIGNAL(actionCheck()), strokeBrush, SIGNAL(actionUncheck()));

	viewHarmonicField = new QAction(tr("&Show harmonic field"), this);
	viewHarmonicField->setIcon(QIcon(":/icons/images/paint.png"));
	viewHarmonicField->setStatusTip(tr("Show harmonic field"));
	connect(viewHarmonicField, SIGNAL(triggered()), m_pViewer, SLOT(showField()));
	connect(viewHarmonicField, SIGNAL(triggered()), this, SLOT(quitSelection()));
	connect(viewHarmonicField, SIGNAL(triggered()), selectAction, SIGNAL(actionUncheck()));

	viewIsolines = new QAction(tr("&Show isolines"), this);
	viewIsolines->setIcon(QIcon(":/icons/images/show-iso.png"));
	viewIsolines->setStatusTip(tr("Show isolines"));
	connect(viewIsolines, SIGNAL(triggered()), m_pViewer, SLOT(showIsolines()));

	viewFinalCut = new QAction(tr("&Show final cut"), this);
	viewFinalCut->setIcon(QIcon(":/icons/images/cut.png"));
	viewFinalCut->setStatusTip(tr("Show final cut"));
	connect(viewFinalCut, SIGNAL(triggered()), m_pViewer, SLOT(showCutResult()));
}

void MainWindow::createToolbar()
{
	fileToolbar = addToolBar(tr("&File"));
	fileToolbar->addAction(openMeshAction);
	fileToolbar->addAction(openTextureAction);

	editToolbar = addToolBar(tr("&Edit"));
	editToolbar->addAction(strokeBrush);
	editToolbar->addAction(selectAction);
	//editToolbar->addWidget(pickColor);
	editToolbar->addAction(viewHarmonicField);
	editToolbar->addAction(viewIsolines);
	editToolbar->addAction(viewFinalCut);
	//editToolbar->addAction(takeSnapshot);

	viewToolbar = addToolBar(tr("&View"));
	viewToolbar->addAction(renderOptions);
	viewToolbar->addAction(viewTexture);
	viewToolbar->addAction(viewPoints);
	viewToolbar->addAction(viewWireframe);
	viewToolbar->addAction(viewFlatlines);
	viewToolbar->addAction(viewFlat);
	viewToolbar->addAction(viewSmooth);
	viewToolbar->addAction(viewBoundary);

	drawModeGroup = new QActionGroup(this);
	drawModeGroup->addAction(viewPoints);
	drawModeGroup->addAction(viewWireframe);
	drawModeGroup->addAction(viewFlatlines);
	drawModeGroup->addAction(viewFlat);
	drawModeGroup->addAction(viewSmooth);
	drawModeGroup->addAction(viewTexture);
}

void MainWindow::setRenderOptions()
{
	RenderOptionsWindow * pWin = new RenderOptionsWindow(m_pViewer,NULL);
	pWin->show();
}

void MainWindow::dragEnterEvent(QDragEnterEvent * e)
{
	if (e->mimeData()->hasUrls())
	{
		e->acceptProposedAction();
	}
}

void MainWindow::dropEvent(QDropEvent * e)
{
	for each (const QUrl &url in e->mimeData()->urls())
	{
		const QString & filename = url.toLocalFile();
		std::string sFileName = filename.toStdString();
		std::string fExt;

		if (filename.endsWith(".m"))
		{
			fExt = "m";
		}
		else if (filename.endsWith(".obj"))
		{
			fExt = "obj";
		}
		else if (filename.endsWith(".bmp"))
		{
			fExt = "bmp";
		}
		else
		{
			QMessageBox dragFileFailed;
			dragFileFailed.setText("Illegal File Type.");
			dragFileFailed.exec();
		}
		m_pViewer->loadFromMainWindow(sFileName, fExt);
	
	}
}

void MainWindow::enterSelection()
{
	selectAction->setChecked(true);
}

void MainWindow::quitSelection()
{
	selectAction->setChecked(false);
}

void MainWindow::enterBrush()
{
	strokeBrush->setChecked(true);
	if (m_pViewer->get_update_mode())
	{
		setCursor(*updateCursor);
	}
	else
	{
		setCursor(Qt::CrossCursor);
	}
}

void MainWindow::quitBrush()
{
	strokeBrush->setChecked(false);
	setCursor(Qt::ArrowCursor);
}

void MainWindow::keyPressEvent(QKeyEvent * keyEvent)
{
	switch (keyEvent->key())
	{
	case Qt::Key_S:
		if (selectAction->isChecked())
		{
			quitSelection();
			selectAction->actionUncheck();
		}
		else
		{
			enterSelection();
			selectAction->actionCheck();
		}
		break;
	case Qt::Key_B:
		if (strokeBrush->isChecked())
		{
			quitBrush();
			strokeBrush->actionUncheck();
		}
		else
		{
			enterBrush();
			strokeBrush->actionCheck();
		}
		break;
	case Qt::Key_U:
		if (strokeBrush->isChecked())
		{
			if (m_pViewer->get_update_mode())
			{
				m_pViewer->set_update_mode(false);
				setCursor(Qt::CrossCursor);
			}
			else
			{
				m_pViewer->set_update_mode(true);
				setCursor(*updateCursor);
			}
		}
		break;
	default:
		break;
	}
}