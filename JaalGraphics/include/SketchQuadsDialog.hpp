#pragma once

#include <QDialog>

#include "Ui_SketchQuadsDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshGeodesics.hpp"
#include "RayTracer.hpp"

class JSketchQuadsDialog : public QDialog, public Ui::SketchQuadsDialog {
    Q_OBJECT

public:
    JSketchQuadsDialog( QWidget *parent = 0);
    ~JSketchQuadsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setMesh( const JMeshPtr &m);

//private slots:
private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

private:
    JMeshPtr mesh;
    boost::scoped_ptr<JRayTracer> rayTracer;
    boost::scoped_ptr<JMeshGeodesics> geoPath;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////

