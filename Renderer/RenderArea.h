//
//  RenderArea.h
//  Renderer
//
//  Created by Raphael Kargon on 6/6/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__RenderArea__
#define __Renderer__RenderArea__

#include <algorithm>
#include <fstream>
#include <QtCore/QTimer>
#include <QtGui/QColor>
#include <QtGui/QKeyEvent>
#include <QtGui/QMouseEvent>
#include <QtGui/QPalette>
#include <QtGui/QPainter>
#include <QtGui/QWheelEvent>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLayout>
#include <QtWidgets/QWidget>
#include "camera.h"
#include "kdtree.h"
#include "mesh.h"
#include "rendering.h"

class RenderArea : public QWidget
{
public:
    std::vector<edge*> kdedges;
    scene *sc;
    
    bool amboc=false;
    QImage *renderimg;
    unsigned int rendermode=0;
    
    real *zbuffer;
    int *normalmap;
    
    QPoint prevpos; //for mouse movement tracking
    QLabel *statuslbl;
    
    static const int tilesize = 32;
    
    RenderArea(QWidget *parent = 0);
    QSize minimumSizeHint() const;
    QSize sizeHint() const;
    
    void updateText();
    void updateImage();
    void drawWireFrame();
    void zBufferDraw();
    
    // The functions below all rely on generate_maps_vector, where the real work is done
    void generate_maps_vector(int mapflags);
    void zBufferDraw_vector();
    void paintNormalMap();
    void SSAO();
    
    void rayTraceUnthreaded();
    
protected:
    void paintEvent(QPaintEvent *event);
    void keyPressEvent(QKeyEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);
    void resizeEvent(QResizeEvent *event);
};

void rayTraceTile(QImage *renderimg);

QColor colorToQColor(const color& c);


//stores increment values of barycentric coordinates for 4 pixel values at once
typedef struct EdgeVect{
    static const int stepX = 4;
    static const int stepY = 1;
    
    __v4si oneStepX;
    __v4si oneStepY;
    
    __v4si init(const point2D<int>& v0, const point2D<int>& v1, const point2D<int>& origin);
} EdgeVect;

#endif /* defined(__Renderer__RenderArea__) */
