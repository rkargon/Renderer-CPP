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
    std::ifstream dragonfile("/Users/raphaelkargon/Documents/Programming/Computer Graphics/Models/dragonsmall.stl");
    std::ifstream spherefile("/Users/raphaelkargon/Documents/Programming/Computer Graphics/Models/sphere.stl");
    camera *cam = new camera();
    cam->dof_focus_distance = 4;
    cam->aperture_size = 0.0;
    std::vector<lamp*> lamps;
//    lamps.push_back(new lamp(45, 2, vertex(-10, 0,-7), rgb_to_color(0xFFAAAA)));
//    lamps.push_back(new lamp(45, 2, vertex( 10, 0,-7), rgb_to_color(0xAAFFAA)));
//    lamps.push_back(new lamp(45, 2, vertex( 0,-10, 7), rgb_to_color(0xAAAAFF)));
    lamps.push_back(new lamp(25, 2, vertex( 0, 5, 5), rgb_to_color(0xFFCC66)));
    world* sc_world = new world(color(0.4,0.4,0.4), color(0.3,0.4,0.5));;
    std::vector<mesh*> objects;

    //dragon
    mesh *dragonobj = new mesh(dragonfile, "Dragon");
    dragonfile.close();
    dragonobj->mat = new material();
    dragonobj->bsdf = new DiffuseBSDF();
    dragonobj->project_texture(TEX_PROJ_SPHERICAL);

    //sphere
    mesh *sphereobj = new mesh(spherefile, "Sphere");
    spherefile.close();
    sphereobj->mat = new material();
    sphereobj->bsdf = new EmissionBSDF();
    sphereobj->project_texture(TEX_PROJ_SPHERICAL);
    sphereobj->scale_centered(vertex(0.3, 0.3, 0.3));
    sphereobj->move(vertex(0, 0, 1));

    objects.push_back(dragonobj);
    objects.push_back(sphereobj);

    distance_estimator *de_obj = new auto(de_mandelbox(2.5, 1, 0.5, 1, 30));
    sc = new scene(cam, lamps, sc_world, objects, de_obj, new material{});
//    if(sc->kdt != nullptr) sc->kdt->print_stats();
    
    //set up images and buffers
    int h = height(), w = width();
    imgrasters = new raster(w, h);
    renderimg = new QImage((uchar*) imgrasters->colbuffer, width(), height(), QImage::Format_ARGB32);
    manager = new thread_manager(&imgrasters, sc, 1);

    //status label
    statuslbl = new QLabel("Raph Renderer 2017");
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
            statuslbl->setText("2. ZBuffer Draw");
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
        case 7:
            statuslbl->setText("7. Fractal Rendering");
            break;
        default:
            statuslbl->setText("Raph Renderer 2015");
    }
}

void RenderArea::updateImage(){
    manager->stop();

    switch(rendermode){
        case 0:
            drawWireFrame();
            break;
        case 1:
            zbuffer_draw(imgrasters, sc);
            break;
        case 2:
            SSAO(imgrasters, sc);
            break;
        case 3:
            paint_normal_map(imgrasters, sc);
            break;
        case 4:
            manager->set_render_method(ray_trace_pixel);
            manager->start();
            break;
        case 5:
            manager->set_render_method(path_trace_pixel);
            manager->start();
            break;
        case 6:
            manager->set_render_method(amb_occ_pixel);
            manager->start();
            break;
        case 7:
            manager->set_render_method(ray_march_pixel);
            manager->start();
            break;
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
            rendermode  = (rendermode+8)%8;
            updateText();
            updateImage();
            break;
        case Qt::Key_Space:
            sc->cam->focus = vertex();
            sc->cam->center_focus();
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
        sc->cam->shift_focus(-delta.x()/100.0, delta.y()/100.0);
    }
    else{
        sc->cam->rotate_local_x(delta.y()/100.0);
        sc->cam->rotate_local_y(delta.x()/100.0);
    }
    sc->cam->center_focus();
    prevpos = pos;
    updateImage();
    repaint();
}

void RenderArea::mousePressEvent(QMouseEvent *event){
    prevpos = event->pos();
    if(event->button()==Qt::MiddleButton){
        sc->cam->set_global_rotation(0,0,0);
        sc->cam->center_focus();
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
    // thread manager will segfault if image buffer changes out from under its feet.
    manager->stop();
    imgrasters->resize(width(), height());
    renderimg = new QImage((uchar*) imgrasters->colbuffer, width(), height(), QImage::Format_RGB32);
    updateImage();
}

/* DRAWING FUNCTIONS */

void RenderArea::drawWireFrame(){
    renderimg->fill(0xffffff);
    QPainter painter(renderimg);
    int w = width(), h = height();
    point_2d<double> p1,p2;
    painter.setPen(QColor(0,0,0));

    for(mesh *obj: sc->objects){
        for(edge *e : obj->edges){
            p1 = sc->cam->project_vertex(*e->v1, w, h);
            p2 = sc->cam->project_vertex(*e->v2, w, h);
            if(isnan(p1.x)||isnan(p2.x)||isnan(p1.y)||isnan(p2.y)) continue;
            painter.drawLine(p1.x, p1.y, p2.x, p2.y);
        }
    }
}

QColor colorToQColor(const color& c){
    return QColor((int)clamp(c.r*255, 0, 255), (int)clamp(c.g*255, 0, 255), (int)clamp(c.b*255, 0, 255));
}
