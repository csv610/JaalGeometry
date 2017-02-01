#pragma once

#include <QDialog>

#include "Ui_MeshNodesDialog.hpp"

#include "MeshViewer.hpp"
#include "PermuteDialog.hpp"
#include "NodeAttributesDialog.hpp"
#include "MeshNormalsDialog.hpp"
#include "MeshSamplesDialog.hpp"
#include "CuthillMcKeeDialog.hpp"

class JMeshEntityAttribListDialog;

class JMeshNodesDialog : public QDialog, public Ui::MeshNodesDialog {
    Q_OBJECT

public:
    JMeshNodesDialog( QWidget *parent = 0);
    void setViewManager( JaalViewer *v)
    {
        viewManager = v;
        init();
    }
    void setMesh(const JMeshPtr &m);


protected:
    void keyPressEvent( QKeyEvent *e);

private slots:
    void checkNodes();
    void getBoundary();
    void setNumVisible();
    void checkDisplay();
    void checkState();

    void setInternal();
    void setBoundary();

    void centerAtNode();
    void geomCenter();
    void closeDialog();
    void setDefaultColorMethod();

    void updateNodes();
    void changeNodeID();
    void showCoords();
    void modifyCoords();
    void setColorScheme();
    void displayIDs();
    void renumber();

    void openSamplesDialog();
    void openAttribListDialog();
    void openPermuteDialog();
    void openNormalsDialog();
    void getRCMOrdering();

    /*
         void alignNormal( bool refresh = 1);
         void changeCenter( bool refresh = 1);
    */
private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    JMeshPtr mesh;
    JNodePtr currNode;
    JNodeDraw    *drawNode;

    boost::scoped_ptr<JNodeAttributesDialog> attribDialog;
    boost::scoped_ptr<JMeshEntityAttribListDialog> attribListDialog;
    boost::scoped_ptr<JPermuteDialog>  permuteDialog;
    boost::scoped_ptr<JMeshNormalsDialog> normalsDialog;
    boost::scoped_ptr<JNodeColor> nodeColor;
    boost::scoped_ptr<JMeshSamplesDialog> samplesDialog;
    boost::scoped_ptr<JCuthillMcKeeDialog> cuthillMcKeeDialog;

    void init();
    void reset_neighs( bool val );
    void activate( const JNodePtr &vtx );
    void makeConnections();
};

