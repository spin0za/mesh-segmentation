#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QAction>
#include <QtGui>
#include <QToolBar>
#include <QtWidgets>
#include <checkableAction.h>
#include "MeshViewer.h"
#include "RenderOptionsWindow.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
	void createActions();
	void createToolbar();
	void dragEnterEvent(QDragEnterEvent * e);
	void dropEvent(QDropEvent * e);

	void keyPressEvent(QKeyEvent * keyEvent);
	//void keyReleaseEvent(QKeyEvent * keyEvent);

private:
	CMeshViewer * m_pViewer;

	QAction * openMeshAction;
	QAction * openTextureAction;

	QAction * renderOptions;
	QAction * viewPoints;
	QAction * viewWireframe;
	QAction * viewFlatlines;
	QAction * viewFlat;
	QAction * viewSmooth;
	QAction * viewBoundary;
	QAction * viewTexture;
	QAction * takeSnapshot;
	QAction * viewHarmonicField;
	QAction * viewIsolines;
	QAction * viewFinalCut;

	QToolBar * fileToolbar;
	QToolBar * viewToolbar;
	QToolBar * editToolbar;

	QActionGroup * drawModeGroup;

	checkableAction * selectAction;
	checkableAction * strokeBrush;
	
	QCursor *updateCursor;

public slots:
	void setRenderOptions();

private slots:
	void enterSelection();
	void quitSelection();
	void enterBrush();
	void quitBrush();
};

#endif // MAINWINDOW_H
