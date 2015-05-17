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
#include <thread>
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
#include "raster.h"
#include "rendering.h"
#include "renderthreads.h"

class RenderArea : public QWidget
{
public:
    scene *sc;
    
    unsigned int rendermode=0;
    raster *imgrasters;
    QImage *renderimg;
    thread_manager *manager;
    
    QPoint prevpos; //for mouse movement tracking
    QLabel *statuslbl;
    
    RenderArea(QWidget *parent = 0);
    //~RenderArea();
    QSize minimumSizeHint() const;
    QSize sizeHint() const;
    
    void updateText();
    void updateImage();
    void drawWireFrame();
    
protected:
    void paintEvent(QPaintEvent *event);
    void keyPressEvent(QKeyEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);
    void resizeEvent(QResizeEvent *event);
};

QColor colorToQColor(const color& c);


#endif /* defined(__Renderer__RenderArea__) */
