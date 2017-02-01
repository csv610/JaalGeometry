#pragma once

#include <QtWidgets>
#include <QGLViewer/qglviewer.h>
#include <qapplication.h>

#include <string>
#include <vector>
#include "JaalViewer.hpp"
#include "Image.hpp"

class JImageViewer : public JViewComponent {
public:
    static JViewComponentPtr registerComponent(JaalViewer *root);

    explicit JImageViewer( JaalViewer *p);

    ~JImageViewer();

    JaalViewer*  getViewManager() const
    {
        return viewManager;
    }

    void draw();
    void refreshDisplay();

    JImagePtr getImage(const string &s) const
    {
        for(const auto &img: imagedb)
            if( img->getName() == s) return img;
        return nullptr;
    }

    JImagePtr getImage(size_t i = 0) const
    {
        if( imagedb.empty() || i >= imagedb.size() ) return nullptr;
        return imagedb[i];
    }

    int  remove(const JImagePtr &img);

    void addImage( const string &str);

    int getSize() const
    {
        return imagedb.size();
    }
private:
    JaalViewer *viewManager;
    vector<JImagePtr> imagedb;
};

typedef boost::shared_ptr<JImageViewer> JImageViewerPtr;
