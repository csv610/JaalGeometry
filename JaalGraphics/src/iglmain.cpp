#include <QtGui>
#include <iostream>
using namespace std;

#include "IGLViewer.hpp"

int main(int argc, char** argv)
{
    QApplication application(argc,argv);

    IGLMeshPtr mesh( new IGLMesh);
    mesh->readFile(argv[1]);

    IGLViewer  viewer;
    viewer.addMesh(mesh);

    viewer.setWindowTitle("IGLViewer");
    viewer.show();
    return application.exec();

}

