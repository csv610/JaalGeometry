//#include <QtGui>
//#include <iostream>
//using namespace std;

#include <QMainWindow>
#include "MainWindow.hpp"


int main(int argc, char** argv)
{
    QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
    QApplication application(argc,argv);
//  QApplication::setStyle(new QPlastiqueStyle);

//  QFont newFont("Comic_Sans_MS", 14, QFont::Light);
//  QFont newFont("helvetica", 15, QFont::Light);

//  QFont newFont("Times", 10);
//    newFont.setStyleStrategy( QFont::PreferAntialias );
//    QApplication::setFont(newFont);

    JaalMainWindow  win;
    win.setWindowTitle("JaalViewer");

    win.show();
    return application.exec();

}

