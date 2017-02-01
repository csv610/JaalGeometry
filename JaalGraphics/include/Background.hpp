#pragma once

#ifndef BACKGROUND_H
#define BACKGROUND_H

class Background
{
public:
    Background();
    Background( Color &c);
    Background( float r, float g, float b );

    Color getColor();
    void setColor( Color &c);
    void setColor( float r, float g, float b );
};

#endif




