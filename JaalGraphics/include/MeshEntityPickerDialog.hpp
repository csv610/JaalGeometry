#pragma once

#include <QDialog>
#include <QColorDialog>
#include <QMessageBox>

#include "MeshViewer.hpp"
#include "Ui_MeshEntityPickerDialog.hpp"

class JMeshEntityPickerDialog : public QDialog, public Ui::MeshEntityPickerDialog
{
    Q_OBJECT

    static const int SELECT_SHOW   = 0;
    static const int SELECT_HIDE   = 1;
    static const int SELECT_DELETE = 2;

public:
    static const int  SINGLE_SELECTION   = 1;
    static const int  MULTIPLE_SELECTION = 2;

    JMeshEntityPickerDialog( QWidget *parent = 0);

    void setViewManager( JaalViewer *v)
    {
        viewManager = v;
        init();
    }

    void setMesh( const JMeshPtr &m);

protected:
    virtual void  mouseReleaseEvent( QMouseEvent *e);

private slots:

    void setPickNodeColor();
    void setPickEdgeColor();
    void setPickFaceColor();
    void setPickCellColor();

    void setPickHeight();
    void setPickWidth();

    void selectEntity();

    void setSelectionMode();
    void clearAll();
    void setMark();

    void displaySelected();
    void setPickStatus();
    void deleteSelected();
    void numPicked();
    void addManually();
    void closeDialog();
    void setSensitivity();
    void setLensRadius();
    void hideSelected();
    void invertSelected();
    void selectAll();
    void growSelected();
    void shrinkSelected();

    void reject();

private:
    JaalViewer *viewManager;

    JMeshViewerPtr meshViewer;
    JMeshPtr mesh;
    boost::shared_ptr<JMeshEntityPicker> entityPicker;
    int   pickable = 0;
    int   selectOp = SELECT_SHOW;

    void init();
    void makeConnections();
};
