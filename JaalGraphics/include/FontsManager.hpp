#pragma once

#include <map>
#include <string>
#include <sstream>
#include <FTGL/ftgl.h>

#include "EntityColor.hpp"

#define DEFAULT_FONT_PATH    "/usr/share/fonts/truetype"
#define DEFAULT_FONT_FAMILY  "ttf-dejavu"
#define DEFAULT_FONT_NAME    "DejaVuSerif.ttf"

//FTTextureFont* myFont = FTGLFontManager::Instance().GetFont("arial.ttf", 72);

using namespace std;

typedef map<string, FTFont*> FontList;
typedef FontList::const_iterator FontIter;

class FontsManager {
public:
    static const int BITMAP_FONT   = 0;
    static const int BUFFER_FONT   = 1;
    static const int EXTRUDE_FONT  = 2;
    static const int OUTLINE_FONT  = 3;
    static const int PIXMAP_FONT   = 4;
    static const int POLYGON_FONT  = 5;
    static const int TEXTURE_FONT  = 6;

    static FontsManager& Instance();
    ~FontsManager();

    void setPath( const string &p) {
        path = p;
    }
    void setFamily( const string &f) {
        family = f;
    }

    FTFont* getFont(const char *filename = NULL, int size = 72, int type =  POLYGON_FONT);

    void setFontScale( double v ) {
        fontScale = v;
    }

    double getFontScale() const {
        return fontScale;
    }

    void setColor( const JColor &c ) {
        color = c;
    }
    const JColor &getColor() const {
        return color;
    }

    void setRotateAngles(Point3F &p) {
        rotAngles = p;
    }
    const Point3F &getRotateAngles() const {
        return rotAngles;
    }

private:
    // Hide these 'cause this is a singleton.
    FontsManager();
    FontsManager(const FontsManager&) {};
    FontsManager& operator = (const FontsManager&) {
        return *this;
    };

    // container for fonts
    FontList fonts;
    string path, family, fontfile;
    float fontScale;
    JColor  color;
    Point3F rotAngles;

};
