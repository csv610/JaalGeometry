#pragma once

#include <QDialog>
#include "JaalViewer.hpp"
#include "Ui_MagnifyingLensDialog.hpp"

class JMagnifyingLensDialog : public QDialog, public Ui::MagnifyingLensDialog {
    Q_OBJECT

public:
    JMagnifyingLensDialog( QWidget *parent = 0);
    ~JMagnifyingLensDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:
    void  updateCamera();
    void  loadCamera();
    void  addNewLens();
    void  captureRegion();
    void  deleteLens();
    void  closeDialog();

protected:
    virtual void mouseMoveEvent(QMouseEvent *e);
    void actionMouseEvent(int id);

private:
    bool locked = 0;
    JaalViewer *viewManager;
    JMagnifyingLensPtr currLens;
    vector<JMagnifyingLensPtr> lensList;
    double borderWidth;
    JColor borderColor;

    void init();
    void makeConnections();
};
