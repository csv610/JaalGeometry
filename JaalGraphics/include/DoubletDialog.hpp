#pragma once

#include <QDialog>
#include <QColorDialog>
#include "Ui_DoubletDialog.hpp"

#include "MeshViewer.hpp"
#include "QuadCleanUp.hpp"

class DoubletColor : public JNodeColor {
public:
    string getName() const {
        return "DoubletColor";
    }

    bool isPerNode() const {
        return 1;
    }

    DoubletColor();

    int assign(const JNodePtr &v);
    int  operator() (const JNodePtr &v) {
        return assign(v);
    }

private:
    JColor nodoubletColor;
};

///////////////////////////////////////////////////////////////////////////////

class JDoubletDialog : public QDialog, public Ui::DoubletDialog {
    Q_OBJECT

public:
    JDoubletDialog( QWidget *parent = 0);
    ~JDoubletDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void setColor();
    void removeAll();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    boost::shared_ptr<JNodeColor> doubletColor;
    JDoublet  doublet;

    void init();
    void makeConnections();
    void assignColor();
};
