#pragma once

#include <QDialog>

#include "Ui_LocallyInjectiveMapParamsDialog.hpp"
#include "LocallyInjectiveMap.hpp"

class JLocallyInjectiveMapParamsDialog : public QDialog, public Ui::LocallyInjectiveMapParamsDialog {
    Q_OBJECT

public:
    JLocallyInjectiveMapParamsDialog( QWidget *parent = 0);
    ~JLocallyInjectiveMapParamsDialog();

    void setDeformer(JLocallyInjectiveMap *d) {
        limDeformer = d;
        init();
    }

private slots:
    void setParams();

private:
    void init();
    void makeConnections();
    JLocallyInjectiveMap *limDeformer;
};

