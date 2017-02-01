#pragma once

#include <QDialog>

#include "Ui_MeshFacesDialog.hpp"

#include "MeshViewer.hpp"

#include "FaceAttributesDialog.hpp"
#include "MeshNormalsDialog.hpp"

//class MeshEntityAttributesDialog;

class JMeshFacesDialog : public QDialog, public Ui::MeshFacesDialog {
    Q_OBJECT

public:
    JMeshFacesDialog( QWidget *parent = 0);

    void setViewManager( JaalViewer *v)
    {
        viewManager = v;
        init();
    }
    void setMesh(const JMeshPtr &m);

private slots:
    void closeDialog();

    void checkDisplay();
    void displayFaces();
    void displayFaceType();

    void setInternal();
    void setBoundary();
    void setInterface();

    void getBoundary();

    void setNumVisible();
    void checkState();

    void keyPressEvent( QKeyEvent *e);

    void deleteMesh();

    void openAttribListDialog();
    void openNormalsDialog();
    void lookAt();
    void checkLights();
    void reverseAll();

//   void setDefaultColorMethod();
//   void alignAlongAxis( bool refresh = 1);
//   void changeCenter( bool refresh = 1);

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    JFaceDraw *drawFace;
    JFacePtr currFace;

    boost::scoped_ptr<JMeshEntityAttribListDialog> attribListDialog;
    boost::scoped_ptr<JFaceColor> faceColor;
    boost::scoped_ptr<JFaceAttributesDialog>  attribDialog;
    boost::scoped_ptr<JMeshNormalsDialog> normalsDialog;

    void init();
    void makeConnections();
};

