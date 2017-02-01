#pragma once

#include <QDialog>

#include "MeshViewer.hpp"

#include "Ui_SuggestiveContoursDialog.hpp"
#include "SuggestiveContoursViewer.hpp"

class JSuggestiveContoursDialog : public QDialog, public Ui::SuggestiveContoursDialog {
    Q_OBJECT

public:
    JSuggestiveContoursDialog( QWidget *parent = 0);
    ~JSuggestiveContoursDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void closeDialog();
    void setParam();
    void smoothNormals();
    void smoothCurvature();
    void smoothCurvatureDeriv();
    void subdivisionSurface();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    boost::shared_ptr<JSuggestiveContoursViewer> scViewer;

private:
    JMeshPtr mesh;
    void init();
    void makeConnections();

};
////////////////////////////////////////////////////////////////////////////////
