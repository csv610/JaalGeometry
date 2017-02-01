#pragma once

#ifndef LINEATTRIB_H
#define LINEATTRIB_H

class LineAttributes
{
public:
    static const int  PATTERN_SOLID         = 0;
    static const int  PATTERN_DASH          = 0;
    static const int  PATTERN_DOT           = 0;
    static const int  PATTERN_DASH_DOT      = 0;
    static const int  PATTERN_USER_DEFINED  = 0;

    LineAttributes();

    void  setWidth( float w );
    float getWidth() const;

    void  setPattern( char *c);
};

#endif
