#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_SingletDialog.hpp"
#include "MeshViewer.hpp"
#include "QuadCleanUp.hpp"

class SingletColor : public JNodeColor {
public:

    string getName() const {
        return "SingletColor";
    }

    bool isPerNode() const {
        return 1;
    }

    SingletColor();

    int  assign(const JNodePtr &v);

    int  operator() (const JNodePtr &v) {
        return assign(v);
    }

private:
    JColor nosingletColor;
};

////////////////////////////////////////////////////////////////////////////////

class JSingletDialog : public QDialog, public Ui::SingletDialog {
    Q_OBJECT

public:
    JSingletDialog( QWidget *parent = 0);
    ~JSingletDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void searchSinglets();
    void setColor();
    void removeAll();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    boost::shared_ptr<JSinglet> singlets;
    boost::shared_ptr<JNodeColor> singletColor;

    void init();
    int  isQuadMesh();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
