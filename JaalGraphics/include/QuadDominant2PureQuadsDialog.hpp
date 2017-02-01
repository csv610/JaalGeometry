#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_QuadDominant2PureQuadsDialog.hpp"
#include "MeshViewer.hpp"
#include "QuadDominant2PureQuadsMesher.hpp"
#include "MeshEntityPicker.hpp"

class JQuadDominant2PureQuadsDialog : public QDialog, public Ui::QuadDominant2PureQuadsDialog {
    Q_OBJECT

public:
    JQuadDominant2PureQuadsDialog( QWidget *parent = 0);
    ~JQuadDominant2PureQuadsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

protected:
    void mouseReleaseEvent(QMouseEvent *e);

private slots:

    void enumFaces();
    void getNewStrip();
    void remeshStrip();

    void getAllStrips();
    void remeshAllStrips();

    void refineAll();
    void clearAll();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    JMeshEntityPickerPtr entityPicker;

    boost::scoped_ptr<JQuadDominant2PureQuadsMesher>  pureQuadsMesher;

    JFaceSequence newStrip;

    JMeshPtr mesh;

    void setColors();
    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
