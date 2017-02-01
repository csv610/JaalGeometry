#pragma once

#include <QDialog>

#include "Ui_NodeAttributesDialog.hpp"
#include "MeshViewer.hpp"

using namespace Jaal;

class JNodeAttributesDialog : public QDialog, public Ui::NodeAttributesDialog {
    Q_OBJECT

public:
    typedef boost::shared_ptr<JNodeAttributesDialog> shared_ptr;

    JNodeAttributesDialog( QWidget *parent = 0);

    void setViewManager(JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh(const JMeshPtr &m);

    void updateDefaults(bool b) {
        update_default_values = b;
    }

    void setNodes(JNodeSequence &seq);

signals:
    void nodesChanged();

protected:
    virtual void keyPressEvent( QKeyEvent *e);

private slots:
    void setColor();
    void setBorderColor();
    void setGlyph();
    void setPointSize();
    void setSphereRadius();
    void setSphResolution();
    void checkDisplay();
    void closeDialog();
    void setBorderThickness();

private:
    JaalViewer  *viewManager;
    JMeshViewerPtr meshViewer;
    JMeshPtr   mesh;

    JNodeSequence  nodes;

    int  node_type;
    bool update_default_values;

    void init();
    void makeConnections();
};
