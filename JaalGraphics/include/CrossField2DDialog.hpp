#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_CrossField2DDialog.hpp"
#include "MeshViewer.hpp"

#include "CrossField2D.hpp"


////////////////////////////////////////////////////////////////////////////////

class JCrossField2DDialog : public QDialog, public Ui::CrossField2DDialog {
    Q_OBJECT

public:
    JCrossField2DDialog( QWidget *parent = 0);
    ~JCrossField2DDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void generateField();
    void alignEdges();
    void drawVectors();
    void showScalarField();
    void setVecLength();
    void setLineWidth();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    JMeshPtr nodeVecField, edgeVecField, faceVecField;
    boost::scoped_ptr<JCrossField2D> crossField;
    double veclen;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
