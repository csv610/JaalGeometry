#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_QuadDiamondsDialog.hpp"

#include "MeshViewer.hpp"
#include "MeshOptDialog.hpp"
#include "MeshViewer.hpp"
#include "QuadCleanUp.hpp"

class JQuadDiamondsDialog : public QDialog, public Ui::QuadDiamondsDialog {
    Q_OBJECT

public:
    JQuadDiamondsDialog( QWidget *parent = 0);
    ~JQuadDiamondsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:
    void setColor();
    void removeAll();
    void meshOpt();
    void incrementalRemove();
    void searchDiamonds();
    void accept();
    void reject();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    boost::shared_ptr<JNodeColor> irregularColor;
    boost::shared_ptr<JFaceColor> diamondColor;

    boost::scoped_ptr<JMeshOptDialog> meshoptDialog;

    void init();
    void makeConnections();
};

