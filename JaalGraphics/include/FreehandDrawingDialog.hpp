#pragma once

#ifndef JFREEHANDDIA_H
#define JFREEHANDDIA_H

#include <QDialog>
#include <QColorDialog>

#include "Ui_FreehandDrawingDialog.hpp"

class JFreehandDrawingDialog : public QDialog, public Ui::FreehandDrawingDialog {
    Q_OBJECT

public:
    JFreehandDrawingDialog( QWidget *parent = 0);
    ~JFreehandDrawingDialog();

private slots:
    void init();
    void accept();
    void reject();

private:
    void makeConnections();
};

#endif




