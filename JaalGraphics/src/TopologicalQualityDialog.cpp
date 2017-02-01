#include "TopologicalQualityDialog.hpp"
#include <iostream>

using namespace std;

///////////////////////////////////////////////////////////////////////////////

JTopologicalQualityDialog :: JTopologicalQualityDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    meshviewer = nullptr;
    currNodeColor = new IrregularNodeColor;
}
///////////////////////////////////////////////////////////////////////////////
JTopologicalQualityDialog :: ~JTopologicalQualityDialog()
{
    delete currNodeColor;
}

void JTopologicalQualityDialog :: makeConnections()
{
}

///////////////////////////////////////////////////////////////////////////////
