#pragma once

#ifndef JPATCHREMESHDIA_H
#define JPATCHREMESHDIA_H

#include <QDialog>
#include <QColorDialog>

#include "MeshViewer.hpp"
#include "QuadCleanUp.hpp"
#include "Ui_PatchQuadmeshingDialog.hpp"

class QuadDefectColor : public JFaceColor {
public:
    string getName() const {
        return "QuadDefectColor";
    }

    bool isPerFace() const {
        return 1;
    }
    void assignColors( const QDefectivePatch *newpatch, const JMeshPtr &mesh );

    int assign( const JFacePtr &face);

    int  operator() (const JFacePtr &f) {
        return assign(f);
    }
};

class JPatchQuadmeshingDialog : public QDialog, public Ui::PatchQuadmeshingDialog {
    Q_OBJECT

public:
    JPatchQuadmeshingDialog( QWidget *parent = 0);

    ~JPatchQuadmeshingDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:

    void remeshAll();
    void setColor();
    void accept();

    void openTemplateDialog();

private:
    JaalViewer *viewManager;
    JMeshViewer *meshViewer;

    JMeshPtr mesh;
    boost::scoped_ptr<JNodeColor>      irregularColor;
    boost::scoped_ptr<QuadDefectColor> patchColor;
    boost::scoped_ptr<QDefectivePatch> newpatch;
    boost::scoped_ptr<JQuadCleanUp>    qClean;
//   boost::scoped_ptr<JQuadmeshingTemplatesDialog> meshingTemplatesDialog;

    void init();
    void makeConnections();
};

#endif
