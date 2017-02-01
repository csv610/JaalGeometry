#pragma once

#include <GL/gl.h>
#include "EntityColor.hpp"

class JMaterial {
public:
    JMaterial() ;
    JMaterial( const string &s) {
        setMaterial(s);
    }

    void setMaterialName();

    void setAmbientColor( Array3F &c, GLenum face = GL_FRONT_AND_BACK);
    void setAmbientColor( Array4F &c, GLenum face = GL_FRONT_AND_BACK);
    JColor getAmbientColor(GLenum face = GL_FRONT_AND_BACK);

    void  setDiffuseColor( Array3F &c, GLenum face = GL_FRONT_AND_BACK);
    void  setDiffuseColor( Array4F &c, GLenum face = GL_FRONT_AND_BACK);
    JColor getDiffuseColor( GLenum face);

    void setSpecularColor( Array3F &c, GLenum face = GL_FRONT_AND_BACK);
    void setSpecularColor( Array4F &c, GLenum face = GL_FRONT_AND_BACK);
    JColor getSpecularColor( GLenum face = GL_FRONT_AND_BACK);

    void setEmissionColor( Array3F &c, GLenum face = GL_FRONT_AND_BACK);
    void setEmissionColor( Array4F &c, GLenum face = GL_FRONT_AND_BACK);
    JColor getEmissionColor( GLenum face = GL_FRONT_AND_BACK);

    void setShininess( float &v, GLenum face = GL_FRONT_AND_BACK);
    float getShininess(GLenum face = GL_FRONT_AND_BACK);

    void setMaterial( const string &n);
    void apply();

private:

    /*
         JMaterial() {
         }

         JMaterial( const JMaterial &) {};
         JMaterial & operator = ( const JMaterial &) {
              return *this;
         }
    */
    void setMaterial( float *v );

    void setDefault();
    void setUserDefined();

    void setEmerald();
    void setJade();
    void setObsidian();
    void setPearl();
    void setRuby();
    void setTurquoise();
    void setBrass();
    void setBronze();
    void setChrome();
    void setCopper();
    void setGold();
    void setSilver();
    void setBlackPlastic();
    void setCyanPlastic();
    void setGreenPlastic();
    void setRedPlastic();
    void setWhitePlastic();
    void setYellowPlastic();
    void setBlackRubber();
    void setCyanRubber();
    void setGreenRubber();
    void setRedRubber();
    void setWhiteRubber();
    void setYellowRubber();

    float  ambient[4], diffuse[4], emission[4], specular[4];
    float  shine;
};
