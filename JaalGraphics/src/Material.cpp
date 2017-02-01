#include "Material.hpp"

JMaterial :: JMaterial()
{
    ambient[0] = 0.2;
    ambient[1] = 0.2;
    ambient[2] = 0.2;
    ambient[3] = 1.0;

    diffuse[0] = 0.8;
    diffuse[1] = 0.8;
    diffuse[2] = 0.8;
    diffuse[3] = 0.8;

    specular[0] = 0.0;
    specular[1] = 0.0;
    specular[2] = 0.0;
    specular[3] = 1.0;

    emission[0] = 0.0;
    emission[1] = 0.0;
    emission[2] = 0.0;
    emission[3] = 1.0;

    shine = 100;

}
void JMaterial:: setAmbientColor( Array3F &c, GLenum face)
{
    ambient[0] = c[0];
    ambient[1] = c[1];
    ambient[2] = c[2];
    glMaterialfv(face, GL_AMBIENT,  ambient);
}

void JMaterial:: setAmbientColor( Array4F &c, GLenum face)
{
    ambient[0] = c[0];
    ambient[1] = c[1];
    ambient[2] = c[2];
    ambient[3] = c[3];
    glMaterialfv(face, GL_AMBIENT,  ambient);
}

void JMaterial:: setDiffuseColor( Array3F &c, GLenum face)
{
    diffuse[0] = c[0];
    diffuse[1] = c[1];
    diffuse[2] = c[2];
    glMaterialfv(face, GL_DIFFUSE,  ambient);
}

void JMaterial:: setDiffuseColor( Array4F &c, GLenum face)
{
    ambient[0] = c[0];
    ambient[1] = c[1];
    ambient[2] = c[2];
    ambient[3] = c[3];
    glMaterialfv(face, GL_DIFFUSE,  ambient);
}

void JMaterial:: setSpecularColor( Array3F &c, GLenum face)
{
    specular[0] = c[0];
    specular[1] = c[1];
    specular[2] = c[2];
    glMaterialfv(face, GL_SPECULAR,  ambient);
}

void JMaterial:: setSpecularColor( Array4F &c, GLenum face)
{
    specular[0] = c[0];
    specular[1] = c[1];
    specular[2] = c[2];
    specular[3] = c[3];
    glMaterialfv(face, GL_SPECULAR,  ambient);
}

void JMaterial:: setEmissionColor( Array3F &c, GLenum face)
{
    emission[0] = c[0];
    emission[1] = c[1];
    emission[2] = c[2];
    glMaterialfv(face, GL_EMISSION,  ambient);
}

void JMaterial:: setEmissionColor( Array4F &c, GLenum face)
{
    emission[0] = c[0];
    emission[1] = c[1];
    emission[2] = c[2];
    emission[3] = c[3];
    glMaterialfv(face, GL_EMISSION,  ambient);
}

/////////////////////////////////////////////////////////////////////////

void JMaterial ::apply()
{
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,  ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,  diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shine );
}

void JMaterial ::setMaterial( float *val )
{
    cout << "HELLO Material " << endl;
    ambient[0] = val[0];
    ambient[1] = val[1];
    ambient[2] = val[2];
    ambient[3] = 1.0;
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,  ambient);

    cout << "HELLO Material " << endl;
    diffuse[0] = val[3];
    diffuse[1] = val[4];
    diffuse[2] = val[5];
    diffuse[3] = 1.0;
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,  diffuse);

    cout << "HELLO Material " << endl;
    specular[0] = val[6];
    specular[1] = val[7];
    specular[2] = val[8];
    specular[3] = 1.0;
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);

    cout << "HELLO Material " << endl;
    shine = val[9];
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shine );
}

void JMaterial ::setMaterial( const string &str )
{
    if( str == "Brass")     setBrass();
    if( str == "Bronze")    setBronze();
    if( str == "Chrome")    setChrome();
    if( str == "Copper")    setCopper();
    if( str == "Emerald")   setEmerald();
    if( str == "Jade")      setJade();
    if( str == "Obsidian")  setObsidian();
    if( str == "Pearl")     setPearl();
    if( str == "Ruby")      setRuby();
    if( str == "Turquoise") setTurquoise();
    if( str == "Default")   setDefault();
    if( str == "UserDefined")  setUserDefined();
}

////////////////////////////////////////////////////////////////////////
void JMaterial :: setUserDefined()
{
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,  ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,  diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shine );
}
////////////////////////////////////////////////////////////////////////

void JMaterial :: setDefault()
{
    static float val[] = {0.2, 0.2, 0.2,
                          0,8, 0.8, 0.8,
                          0.0, 0.0, 0.0,
                          0.0
                         };
    setMaterial(val);
}
////////////////////////////////////////////////////////////////////////

void JMaterial :: setEmerald()
{

    static float val[] = { 0.02150, 0.17450, 0.0215,
                           0.07568, 0.61424, 0.07568,
                           0.633,   0.727811, 0.633,
                           0.6*128
                         };
    setMaterial(val);
}

////////////////////////////////////////////////////////////////////////


void JMaterial :: setJade()
{
    static float val[] = { 0.135, 0.2225, 0.1575,
                           0.54, 0.89,	0.63,
                           0.316228, 0.316228, 0.316228,
                           0.1*128
                         };
    setMaterial(val);
}

void JMaterial :: setObsidian()
{
    static float val[] = { 0.05375, 0.05, 0.06625,
                           0.18275, 0.17, 0.22525,
                           0.332741, 0.328634,	0.346435,
                           0.3*128
                         };
    setMaterial(val);
}

void JMaterial :: setPearl()
{
    static float val[] = { 0.25, 0.20725, 0.20725,
                           1.0, 0.829, 0.829,
                           0.296648, 0.296648,	0.296648,
                           0.088*128
                         };
    setMaterial(val);
}

void  JMaterial :: setRuby()
{

    static float val[] = { 0.1745, 0.01175, 0.01175,
                           0.61424, 0.04136, 0.04136,
                           0.727811,  0.626959, 0.626959,
                           0.6*128
                         };
    setMaterial(val);

}

void JMaterial :: setTurquoise()
{
    static float val[] = { 0.1, 0.18725, 0.1745,
                           0.396, 0.74151, 0.69102,
                           0.297254, 0.30829, 0.306678,
                           0.1*128
                         };
    setMaterial(val);
}

void JMaterial :: setBrass()
{
    static float val[] =  { 0.329412, 0.223529, 0.027451,
                            0.780392, 0.568627, 0.113725,
                            0.992157, 0.941176, 0.807843,
                            0.21794872*128
                          };
    setMaterial(val);
}


void JMaterial :: setBronze()
{
    static float val[] = { 0.2125, 0.1275, 0.054,
                           0.714, 0.4284,0.18144,
                           0.393548, 0.271906,	0.166721,
                           0.2*128
                         };
    setMaterial(val);
}

void JMaterial :: setChrome ()
{
    static float val[] = { 0.25, 0.25,	0.25,
                           0.4, 0.4, 0.4,
                           0.774597, 0.774597, 0.774597,
                           0.6*128
                         };
    setMaterial(val);
}


void JMaterial :: setCopper()
{
    static float val[] = {  0.19125, 0.0735,	0.0225,
                            0.7038, 0.27048, 0.0828,
                            0.256777, 0.137622, 0.086014,
                            0.1*128
                         };
    setMaterial(val);
}

void JMaterial :: setGold()
{
    static float val[] = { 0.24725, 0.1995, 0.0745,
                           0.75164, 0.60648, 0.22648,
                           0.628281, 0.555802, 0.366065,
                           0.4*128
                         };

    setMaterial(val);
}

void JMaterial :: setSilver()
{
    static float val[] = { 0.19225, 0.19225, 	0.19225,
                           0.50754, 0.50754, 	0.50754,
                           0.508273, 0.508273,	0.508273,
                           0.4*128
                         };
    setMaterial(val);
}

void JMaterial :: setBlackPlastic()
{
    static float val[] = {0.0, 0.0, 0.0,
                          0.01, 0.01, 0.01,
                          0.50, 0.50, 0.50,
                          0.25*128
                         };
    setMaterial(val);
}

void JMaterial :: setCyanPlastic()
{
    static float val[] = { 0.0, 0.1, 0.06,
                           0.0, 0.50980392, 0.50980392,
                           0.50196078,	0.50196078, 0.50196078,
                           0.25*128
                         };
    setMaterial(val);
}


void JMaterial :: setGreenPlastic()
{
    static float val[] = { 0.0, 0.0, 0.0,
                           0.1, 0.35, 0.1,
                           0.45, 0.55, 0.45,
                           0.25*128
                         };
    setMaterial(val);
}


void JMaterial :: setRedPlastic()
{
    static float val[] = { 0.0, 0.0, 0.0,
                           0.5, 0.0, 0.0,
                           0.7, 0.6, 0.6,
                           0.25*128
                         };
    setMaterial(val);
}

void JMaterial :: setWhitePlastic()
{
    static float val[] = { 0.0, 0.0, 0.0,
                           0.55, 0.55, 0.55,
                           0.70, 0.70, 0.70,
                           0.25*128
                         };
    setMaterial(val);
}

void JMaterial :: setYellowPlastic()
{
    static float val[] = { 0.0, 0.0, 0.0,
                           0.5, 0.5, 0.0,
                           0.60, 0.60, 0.50,
                           0.25*128
                         };
    setMaterial(val);
}

void JMaterial :: setBlackRubber()
{
    static float val[] = { 0.02, 0.02, 0.02,
                           0.01, 0.01, 0.01,
                           0.4,  0.4, 0.4,
                           0.078125*128
                         };
    setMaterial(val);
}

void JMaterial :: setCyanRubber()
{
    static float val[] = { 0.0, 0.05, 0.05,
                           0.4, 0.5, 0.5,
                           0.04,0.7, 0.7,
                           0.078125*128
                         };
    setMaterial(val);
}

void JMaterial :: setGreenRubber()
{
    static float val[] = { 0.0, 0.05, 0.0,
                           0.4, 0.5, 0.4,
                           0.04, 0.7, 0.04,
                           0.078125*128
                         };
    setMaterial(val);
}

void JMaterial :: setRedRubber()
{
    static float val[] = { 0.05, 0.0, 0.0,
                           0.5, 0.4, 0.4,
                           0.7, 0.04, 0.04,
                           0.078125*128
                         };
    setMaterial(val);
}

void JMaterial :: setWhiteRubber()
{
    static float val[] = { 0.05, 0.05, 0.05,
                           0.5, 0.5, 0.5,
                           0.7, 0.7, 0.7,
                           0.078125*128
                         };
    setMaterial(val);
}

void JMaterial :: setYellowRubber()
{
    static float val[] = { 0.05, 0.05, 0.0,
                           0.5, 0.5, 0.4,
                           0.7, 0.7, 0.04,
                           0.078125*128
                         };
    setMaterial(val);
}

