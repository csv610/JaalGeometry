#pragma once

#include <GL/gl.h>

// Note: Copied from Java3D API


class TransparencyAttributes
{
public:

    // Transparency mode - defines how transparency is applied to this
    // Appearance component object:
    static const int  OPAQUE      = 0;
    static const int  FASTEST     = 1;
    static const int  NICEST      = 2;
    static const int  SCREENDOOR  = 3;
    static const int  BLENDED     = 4;

    // Blend function - used in blended transparency and antialiasing operations.
    // The source function specifies the factor that is multiplied by the source
    // color. This value is added to the product of the destination factor and the
    // destination color. The default source blend function is BLEND_SRC_ALPHA.
    // The source blend function is one of the following:
    static const int  BLEND_ZERO      = 0;
    static const int  BLEND_ONE       = 1;
    static const int  BLEND_SRC_ALPHA = 2;
    static const int  BLEND_ONE_MINUS_SRC_ALPHA = 3;

    // Constructs a TransparencyAttributes object with default parameters. The
    // default values are as follows:
    // transparency mode : NONE
    // transparency value : 0.0
    // source blend function : BLEND_SRC_ALPHA
    // destination blend function : BLEND_ONE_MINUS_SRC_ALPHA

    TransparencyAttribute();
    TransparencyAttribute(int tmode, float val);
    TransparencyAttribute(int tmode, float val,
                          int srcBlendFunc, int dstBlendFunc);

    void  setMode( int m );
    int   getMode() const;

    void  setDstBlendFunction( int df);
    int   getDstBlendFunction() const;

    void  setSrcBlendFunction( int sf);
    int   getSrcBlendFunction() const;

    // Blend value - the amount of transparency to be applied to this Appearance
    // component object. The transparency values are in the range [0.0, 1.0], with
    // 0.0 being fully opaque and 1.0 being fully transparent.

    void  setTransparency( float a);
    float getTransparency() const;
private:
    int   mode;
    int   srcfunc, dstfunc;
    float alpha;
};
