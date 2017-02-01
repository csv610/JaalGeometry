#pragma once

#include <QDialog>

#include "Ui_EdgeAttributesDialog.hpp"
#include "MeshViewer.hpp"

using namespace Jaal;

class JEdgeAttributesDialog : public QDialog, public Ui::EdgeAttributesDialog {
    Q_OBJECT

public:
    JEdgeAttributesDialog( QWidget *parent = 0);

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh(const JMeshPtr &m);

    void updateDefaultValues( bool b ) {
        update_default_values = b;
    }

    void setEdges( const JEdgeSequence &es);

protected:
    virtual void  mousePressEvent( QMouseEvent *e);
    virtual void  mouseReleaseEvent( QMouseEvent *e);

private slots:
    void setColor();
    void setAlpha();
    void setGlyph();
    void setLineWidth();
    void setCylinderRadius();
    void setCylinderSides();
    void checkDisplay();

    void keyPressEvent( QKeyEvent *e);
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    JMeshEntityPickerPtr picker;
    JMeshPtr   mesh;

    JEdgeSequence edges;
    bool update_default_values;

    void init();
    void makeConnections();

};
