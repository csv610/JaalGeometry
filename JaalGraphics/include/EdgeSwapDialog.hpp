#pragma once

#include <QDialog>
#include <QColorDialog>

#include "MeshViewer.hpp"
#include "QuadCleanUp.hpp"
#include "SwapEdges.hpp"

#include "Ui_EdgeSwapDialog.hpp"

class JEdgeSwapDialog : public QDialog, public Ui::EdgeSwapDialog {
    Q_OBJECT

public:
    JEdgeSwapDialog( QWidget *parent = 0);
    ~JEdgeSwapDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:

    void startswap();
    void accept();
    void selectiveEdge();
    void interactiveEdge();
    void swapSelectiveEdge();
    void swapCurrEdge();
    void skipCurrEdge();

private:
    JaalViewer  *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshEntityPicker *entityPicker;

    JMeshPtr mesh;
    JEdgePtr curredge;
    size_t currPos;
    boost::scoped_ptr<JNodeColor>    irregularColor;
    boost::scoped_ptr<JSwapQuadEdge> qedgeswapper;

    void init();
    void makeConnections();
};
