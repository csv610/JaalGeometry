#include "ScreenShotDialog.hpp"


///////////////////////////////////////////////////////////////////////////////

JScreenShotDialog :: JScreenShotDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    viewManager = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JScreenShotDialog :: ~JScreenShotDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JScreenShotDialog :: setViewManager( JaalViewer *v)
{
    viewManager = v;
    prevWidth  =  viewManager->width();
    prevHeight = viewManager->height();
    widthLineEdit->setText( QString::number(prevWidth));
    heightLineEdit->setText( QString::number(prevHeight));
}
///////////////////////////////////////////////////////////////////////////////

void JScreenShotDialog :: getShot()
{
    QString qstr;
    qstr = widthLineEdit->text();
    int w = qstr.toInt();

    qstr  = heightLineEdit->text();
    int h = qstr.toInt();
    viewManager->resize(w,h);
    viewManager->refreshDisplay();

    static QString lastSelectedDirectory;
    qstr  = QFileDialog::getSaveFileName(this,
                                         *new QString("Select File Name "),
                                         lastSelectedDirectory,
                                         *new QString( "Image Format (*.png, .eps, *.pdf, *.ps)"));

    string  file = StdString(qstr);

    if( file.find(".pdf") != string::npos) {
        int state = GL2PS_OVERFLOW, buffsize = 0;
        GLint viewport[4];
        viewport[0] = 0;
        viewport[1] = 0;
        viewport[2] = viewManager->width();
        viewport[3] = viewManager->height();

        FILE *fp = fopen(file.c_str(), "wb");

        if(!fp) {
            cout << "Unable to open file %s for writing" << file << endl;
            return;
        }

        int format  = GL2PS_PDF;
//      int sort    = GL2PS_SIMPLE_SORT;
        int sort    = GL2PS_BSP_SORT;
//      int options = GL2PS_OCCLUSION_CULL | GL2PS_SIMPLE_LINE_OFFSET | GL2PS_SILENT | GL2PS_NO_BLENDING;
        int options = GL2PS_OCCLUSION_CULL | GL2PS_SIMPLE_LINE_OFFSET | GL2PS_SILENT;
        int nbcol   = 0;

        JWaitCursor waitCursor;
        waitCursor.start();
        while(state == GL2PS_OVERFLOW) {
            buffsize += 1024*1024;
            gl2psBeginPage(file.c_str(), "ChamanSinghVerma", viewport, format, sort, options,
                           GL_RGBA, 0, NULL, nbcol, nbcol, nbcol,
                           buffsize, fp, NULL);
            viewManager->refreshDisplay();
            state = gl2psEndPage();
        }
    }

    if( file.find(".png") != string::npos) {
        viewManager->setSnapshotQuality(100);
        viewManager->setSnapshotFormat("PNG");
        viewManager->saveSnapshot( qstr );
    }

    if( file.find(".eps") != string::npos) {
        viewManager->setSnapshotQuality(90);
        viewManager->setSnapshotFormat("EPS");
        viewManager->saveSnapshot( qstr );
    }

    if( file.find(".ps") != string::npos) {
        viewManager->setSnapshotQuality(100);
        viewManager->setSnapshotFormat("PS");
        viewManager->saveSnapshot( qstr);
    }

    viewManager->resize(prevWidth,prevHeight);
    viewManager->updateGL();
}
///////////////////////////////////////////////////////////////////////////////

void JScreenShotDialog :: makeConnections()
{
    connect( getShotPushButton, &QPushButton::clicked,  [=] {getShot();});
    connect( closePushButton, &QPushButton::clicked,  [=] {close();});
}
