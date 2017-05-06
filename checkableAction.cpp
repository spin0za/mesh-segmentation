#include "checkableAction.h"

checkableAction::checkableAction(QObject *parent) : QAction(parent)
{
	connect(this, SIGNAL(triggered()), this, SLOT(forwardCheckSignal()));
}


checkableAction::~checkableAction()
{
}

void checkableAction::forwardCheckSignal()
{
	if (! isCheckable())
	{
		triggered();
	}
	else
	{
		if (isChecked())
		{
			emit actionCheck();
		}
		else
		{
			emit actionUncheck();
		}
	}
}
