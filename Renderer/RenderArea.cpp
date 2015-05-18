//
//  RenderArea.cpp
//  Renderer
//
//  Created by Raphael Kargon on 6/6/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "RenderArea.h"

RenderArea::RenderArea(QWidget *parent): QWidget(parent){
    setBackgroundRole(QPalette::Base);
    setAutoFillBackground(true);
    setFocusPolicy(Qt::StrongFocus);
    setFocus(Qt::ActiveWindowFocusReason);
    setLayout(new QHBoxLayout);
    
    //scene setup
    std::ifstream dragonfile("/Users/raphaelkargon/Documents/Programming/STL Renderer/dragonsmall.stl");
    std::ifstream spherefile("/Users/raphaelkargon/Documents/Programming/STL Renderer/sphere.stl");
    camera *cam = new camera();
    std::vector<lamp*> lamps;
    lamps.push_back(new lamp(15, 2, vertex(-4, 0,-2.828), RGBToColor(0xFFAAAA)));
    lamps.push_back(new lamp(15, 2, vertex( 4, 0,-2.828), RGBToColor(0xAAFFAA)));
    lamps.push_back(new lamp(15, 2, vertex( 0,-4, 2.828), RGBToColor(0xAAAAFF)));
    lamps.push_back(new lamp(15, 2, vertex( 0, 4, 2.828), RGBToColor(0xFFFFAA)));
    world* sc_world = new sky();
    std::vector<mesh*> objects;
    
    //dragon
    mesh *dragonobj = new mesh(dragonfile, "Dragon");
    dragonfile.close();
    dragonobj->mat = new material();
    dragonobj->bsdf = new GlossyBSDF();
    dragonobj->project_texture(TEX_PROJ_SPHERICAL);

    //sphere
    mesh *sphereobj = new mesh(spherefile, "Sphere");
    spherefile.close();
    sphereobj->mat = new material();
    sphereobj->bsdf = new EmissionBSDF(vertex(1,1,1), 15);
    sphereobj->project_texture(TEX_PROJ_SPHERICAL);
    sphereobj->scale_centered(vertex(0.5, 0.5, 0.5));
    sphereobj->move(vertex(5, 0, 0));
    
    objects.push_back(dragonobj);
    objects.push_back(sphereobj);
    sc = new scene(cam, lamps, sc_world, objects);
    if(sc->kdt != nullptr) sc->kdt->printstats();
    
    //set up images and buffers
    int h = height(), w = width();
    imgrasters = new raster(w, h);
    renderimg = new QImage((uchar*) imgrasters->colbuffer, width(), height(), QImage::Format_ARGB32);
    manager = new thread_manager(6, &imgrasters, sc);
    
    //status label
    statuslbl = new QLabel("Raph Renderer 2015");
    statuslbl->setParent(this);
    this->layout()->addWidget(statuslbl);
    statuslbl->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
    statuslbl->setAlignment(Qt::AlignTop);
    updateText();
    statuslbl->show();
}

//RenderArea::~RenderArea(){
//    //TODO actually set up destructors for 3d data
//    delete sc;
//    delete imgrasters;
//    delete renderimg;
//    delete manager;
//}

QSize RenderArea::minimumSizeHint() const{
    return QSize(400,400);
}

QSize RenderArea::sizeHint() const{
    return QSize(10000,10000);
}

void RenderArea::paintEvent(QPaintEvent *event){
    QPainter painter(this);
    painter.drawImage(QPoint(0,0), *renderimg);
}

void RenderArea::updateText(){
    switch(rendermode){
        case 0:
            statuslbl->setText("1. Wireframe");
            break;
        case 1:
            statuslbl->setText("2. ZBuffer draw, SSE");
            break;
        case 2:
            statuslbl->setText("3. SSAO");
            break;
        case 3:
            statuslbl->setText("4. Normal map");
            break;
        case 4:
            statuslbl->setText("5. Raytracing");
            break;
        case 5:
            statuslbl->setText("6. Path Tracing");
            break;
        case 6:
            statuslbl->setText("6. Ambient Occlusion");
            break;
        default:
            statuslbl->setText("Raph Renderer 2015");
    }
}

void RenderArea::updateImage(){
    manager->stop();
    manager->fill_tile_queue();
    
    switch(rendermode){
        case 0:
            drawWireFrame();
            break;
        case 1:
            zBufferDraw_vector(imgrasters, sc);
            break;
        case 2:
            SSAO(imgrasters, sc);
            break;
        case 3:
            paintNormalMap(imgrasters, sc);
            break;
        case 4:
            manager->set_render_method(rayTracePixel);
            manager->start();
            break;
        case 5:
            manager->set_render_method(pathTracePixel);
            manager->start();
            break;
        case 6:
            manager->set_render_method(ambOccPixel);
            manager->start();
            break;
    }
}

/* INPUT LISTENERS */

void RenderArea::keyPressEvent(QKeyEvent *event){
    switch(event->key()){
        case Qt::Key_5:
            sc->cam->ortho = !sc->cam->ortho;
            updateImage();
            break;
        case Qt::Key_S:
            //TODO: add per-object smoothing and, in general, selection
            for(mesh *o : sc->objects){
                o->smooth = !o->smooth;
            }
            updateText();
            updateImage();
            break;
        case Qt::Key_Z:
            if(event->modifiers() & Qt::ShiftModifier) --rendermode;
            else ++rendermode;
            rendermode  = (rendermode+7)%7;
            updateText();
            updateImage();
            break;
        case Qt::Key_Space:
            sc->cam->focus = vertex();
            sc->cam->centerFocus();
            updateImage();
            break;
        default:
            QWidget::keyPressEvent(event);
    }
    repaint();
}

void RenderArea::mouseMoveEvent(QMouseEvent *event){
    QPoint pos = event->pos();
    QPoint delta = pos - prevpos; //prevpos initialized in mousePressEvent
    
    if(event->modifiers() & Qt::ShiftModifier){
        sc->cam->shiftFocus(-delta.x()/100.0, delta.y()/100.0);
    }
    else{
        sc->cam->rotateLocalX(delta.y()/100.0);
        sc->cam->rotateLocalY(delta.x()/100.0);
    }
    sc->cam->centerFocus();
    prevpos = pos;
    updateImage();
    repaint();
}

void RenderArea::mousePressEvent(QMouseEvent *event){
    prevpos = event->pos();
    if(event->button()==Qt::MiddleButton){
        sc->cam->setGlobalRotation(0,0,0);
        sc->cam->centerFocus();
        updateImage();
        repaint();
    }
}

void RenderArea::wheelEvent(QWheelEvent *event){
    double zoomfactor = pow(1.005, -event->delta());
    sc->cam->zoom((double)zoomfactor);
    updateImage();
    repaint();
}

void RenderArea::resizeEvent(QResizeEvent *event){
    imgrasters->resize(width(), height());
    renderimg = new QImage((uchar*) imgrasters->colbuffer, width(), height(), QImage::Format_RGB32);
    updateImage();
}

/* DRAWING FUNCTIONS */

void RenderArea::drawWireFrame(){
    renderimg->fill(0xffffff);
    QPainter painter(renderimg);
    int w = width(), h = height();
    point2D<double> p1,p2;
    painter.setPen(QColor(0,0,0));
    
    for(mesh *obj: sc->objects){
        for(edge *e : obj->edges){
            p1 = sc->cam->projectVertex(*e->v1, w, h);
            p2 = sc->cam->projectVertex(*e->v2, w, h);
            if(isnan(p1.x)||isnan(p2.x)||isnan(p1.y)||isnan(p2.y)) continue;
            painter.drawLine(p1.x, p1.y, p2.x, p2.y);
        }
    }
}

QColor colorToQColor(const color& c){
    return QColor((int)clamp(c.r*255, 0, 255), (int)clamp(c.g*255, 0, 255), (int)clamp(c.b*255, 0, 255));
}