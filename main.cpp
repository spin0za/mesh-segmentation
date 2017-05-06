#include "MainWindow.h"
#include <QApplication>
#include <string>

int main(int argc, char *argv[])
{
	QApplication Viewer(argc, argv);
	MainWindow mainWindow;

	//mainWindow.resize(1024, 768);
	QRect available_geom = QDesktopWidget().availableGeometry();
	int W = available_geom.width();
	int H = available_geom.height();
	int w = 640;
	int h = 480;
	int x = (W - w) / 2;
	int y = (H - h) / 2;
	mainWindow.setWindowIcon(QIcon("david.ico"));
	mainWindow.setGeometry(x, y, w, h);
	mainWindow.showNormal();

	std::string title = "Mesh Cutter Console";
	SetConsoleTitleA(title.c_str());

	return Viewer.exec();
}