#include "FontsManager.hpp"

FontsManager& FontsManager :: Instance()
{
    static FontsManager tm;
    return tm;
}

FontsManager :: FontsManager()
{
    path     = DEFAULT_FONT_PATH;
    family   = DEFAULT_FONT_FAMILY;
    fontfile = DEFAULT_FONT_NAME;

    fontScale = 1.0;

    color[0]  = 1.0;
    color[1]  = 0.5;
    color[2]  = 0.5;
    color[3]  = 1.0;
    rotAngles[0] = 0.0;
    rotAngles[1] = 0.0;
    rotAngles[2] = 0.0;
}

FontsManager ::  ~FontsManager()
{
    FontIter font;
    for(font = fonts.begin(); font != fonts.end(); font++) {
        delete (*font).second;
    }
    fonts.clear();
}

FTFont* FontsManager :: getFont(const char *filename, int size, int type)
{
    ostringstream oss;
    if( filename != nullptr ) fontfile = string(filename);

    oss << path <<  "/" << family << "/" << fontfile << size;
    string fontKey = oss.str();

    switch(type) {
    case BITMAP_FONT:
        fontKey += "Bitmap";
        break;
    case BUFFER_FONT:
        fontKey += "Buffer";
        break;
    case EXTRUDE_FONT:
        fontKey += "Extrude";
        break;
    case OUTLINE_FONT:
        fontKey += "Online";
        break;
    case PIXMAP_FONT:
        fontKey += "Pixmap";
        break;
    case POLYGON_FONT:
        fontKey += "Polygon";
        break;
    case TEXTURE_FONT:
        fontKey += "Texture";
        break;
    }

    FontIter result = fonts.find(fontKey);
    if(result != fonts.end()) {
        return result->second;
    }

    string fullname = path + "/" + family + "/" + fontfile;

    FTFont* font = nullptr;
    switch(type) {
    case BITMAP_FONT:
        font = new FTBitmapFont(fullname.c_str() );
        break;
    case BUFFER_FONT:
        font = new FTBufferFont(fullname.c_str() );
        break;
    case EXTRUDE_FONT:
        font = new FTExtrudeFont(fullname.c_str());
        break;
    case OUTLINE_FONT:
        font = new FTOutlineFont(fullname.c_str());
        break;
    case PIXMAP_FONT:
        font = new FTPixmapFont(fullname.c_str());
        break;
    case POLYGON_FONT:
        font = new FTPolygonFont(fullname.c_str());
        break;
    case TEXTURE_FONT:
        font = new FTTextureFont(fullname.c_str());
        break;
    }

    if( font == nullptr ) return nullptr;


    if( font->Error() ) {
        delete font;
        return nullptr;
    }

    if(!font->FaceSize(size)) {
        delete font;
        return nullptr;
    }

    font->UseDisplayList(1);
    font->CharMap(ft_encoding_unicode);

    fonts[fontKey] = font;
    return font;
}

