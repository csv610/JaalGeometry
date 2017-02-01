#pragma once

#include <QDialog>

#include "Ui_CongruentTessellationDialog.hpp"
#include "MeshViewer.hpp"

///////////////////////////////////////////////////////////////////////////////

class JCongruentViewer : public JViewComponent
{
public :
    JCongruentViewer();

    void setParamCorners( const JNodeSequence &n) {
        paramCorners = n;
    }
    void setFirstNode(JNodePtr f) {
        firstNode = f;
    }
    void setArrows( const JNodeSequence &n) {
        srcdstNodes = n;
    }

    void draw();

private:
    JMeshPtr  mesh;
    JNodeSequence paramCorners;
    JNodePtr firstNode;
    JNodeSequence srcdstNodes;
};

///////////////////////////////////////////////////////////////////////////////


class JCongruentTessellationDialog : public QDialog, public Ui::CongruentTessellationDialog
{
    Q_OBJECT

public:
    JCongruentTessellationDialog( QWidget *parent = 0);
    ~JCongruentTessellationDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

protected:
    virtual void  mouseReleaseEvent( QMouseEvent *e);

private slots:
    void generate();
    void setParamShape();
    void closeDialog();

private:
    JaalViewer *viewManager;

    JMeshViewerPtr meshViewer;
    JEdgeSequence boundedges;
    JNodeSequence paramCorners;
    JNodeSequence srcdstNodes;
    JMeshEntityPickerPtr picker;

    boost::shared_ptr<JCongruentViewer> congruentViewer;

    JMeshPtr mesh;
    void setBoundingSquare();
    void setBoundingTriangle();

    void init();
    void makeConnections();

};





