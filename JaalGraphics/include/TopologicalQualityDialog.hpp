#ifndef JTOPOQUAL_H
#define JTOPOQUAL_H

fdgfdgdf

#include "Ui_TopologicalQualityDialog.hpp"

#include <QDialog>
#include <QColorDialog>
#include "MeshViewer.hpp"


class JTopologicalQualityDialog : public QDialog, public Ui::TopologicalQualityDialog
{
    Q_OBJECT

public:
    JTopologicalQualityDialog( QWidget *parent = 0);
    ~JTopologicalQualityDialog();

    void setMeshViewer( JMeshViewer *v) {
        meshviewer = v;
        meshviewer->setNodeColorMethod( currNodeColor );
        Mesh *mesh = meshviewer->getMesh();
        mesh->build_relations(0,2);
    }

// private slots:

private:
    JMeshViewer *meshviewer;
    NodeColor   *prevNodeColor, *currNodeColor;

    void makeConnections();
};

#endif




