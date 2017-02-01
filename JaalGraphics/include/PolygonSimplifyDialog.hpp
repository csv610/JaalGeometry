#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_PolygonSimplifyDialog.hpp"
#include "MeshViewer.hpp"
//#include "PolygonSimplify.hpp"

class JPolygonSimplifyDialog : public QDialog, public Ui::PolygonSimplifyDialog {
    Q_OBJECT

public:
    JPolygonSimplifyDialog( QWidget *parent = 0);
    ~JPolygonSimplifyDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void getSimplified();
    void acceptMesh();
    void displayMesh();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh, simplifiedMesh;
 // boost::scoped_ptr<JPolygonSimplify> simplify;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
