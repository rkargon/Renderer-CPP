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
    
    mesh *dragonobj = new mesh(dragonfile, "Dragon");
    dragonobj->mat = new material();
    dragonobj->bsdf = new DiffuseBSDF();
    dragonobj->project_texture(TEX_PROJ_SPHERICAL);
    
    mesh *sphereobj = new mesh(spherefile, "Sphere");
    sphereobj->mat = new material();
    //sphereobj->bsdf = new DiffuseBSDF();
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
    
    //status label
    statuslbl = new QLabel("Raph Renderer 2015");
    statuslbl->setParent(this);
    this->layout()->addWidget(statuslbl);
    statuslbl->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
    statuslbl->setAlignment(Qt::AlignTop);
    updateText();
    statuslbl->show();
}

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
            statuslbl->setText(QString("5. Raytracing")  + QString(amboc ? ", ambient occlusion" : ""));
            break;
        case 5:
            statuslbl->setText("6. Path Tracing");
            break;
        case 6:
            statuslbl->setText("7. Stereogram");
            break;
        default:
            statuslbl->setText("Raph Renderer 2015 Happy New Year!");
    }
}

void RenderArea::updateImage(){
    switch(rendermode){
        case 0:
            drawWireFrame();
            break;
        case 1:
            zBufferDraw_vector();
            break;
        case 2:
            SSAO();
            break;
        case 3:
            paintNormalMap();
            break;
        case 4:
            rayTraceUnthreaded();
            break;
        case 5:
            pathTraceUnthreaded();
            break;
        case 6:
            drawStereoGram();
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
        case Qt::Key_O:
            amboc = !amboc;
            updateText();
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
    painter.setPen(QColor(255,0,0));
    for(edge *e : kdedges){
        p1 = sc->cam->projectVertex(*e->v1, w, h);
        p2 = sc->cam->projectVertex(*e->v2, w, h);
        if(isnan(p1.x)||isnan(p2.x)||isnan(p1.y)||isnan(p2.y)) continue;
        painter.drawLine(p1.x, p1.y, p2.x, p2.y);
    }
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

void RenderArea::zBufferDraw_vector(){
    generate_maps_vector(5, imgrasters, sc);
}

void RenderArea::paintNormalMap(){
    generate_maps_vector(3, imgrasters, sc);
    std::fill_n(imgrasters->colbuffer, width()*height(), 0xffffff);
    std::copy(imgrasters->normbuffer, imgrasters->normbuffer+(width()*height()), imgrasters->colbuffer);
}

//not really SSAO, should probably fix some tweaks
void RenderArea::SSAO()
{
    int w=width(), h=height();
    generate_maps_vector(7, imgrasters, sc); //generate depth and normal maps
    //renderimg->fill(0xffffff);
    double ao, d, z, ztmp;
    vertex v, dv, vtmp, n;
    ray r, rtmp;
    int x, y, dy, dx, dir, nsamples;
    int dx_arr[4] = {-3, 3, 0, 0};
    int dy_arr[4] = {0, 0, -3, 3};
    
    for(y=0; y<h; y++){
        for(x=0; x<w; x++){
            z = imgrasters->zbuffer[y*w + x];
            if(z==1) continue;
            z = z * (sc->cam->maxdist-sc->cam->mindist) + sc->cam->mindist;
            r = sc->cam->castRay(x, y, w, h);
            v = r.org + r.dir*z;
            n = RGBToNormal(imgrasters->normbuffer[y*w+x]);
            
            ao = 0;
            nsamples=0;
            for(dir=0; dir<4; dir++){
                dx = dx_arr[dir];
                dy = dy_arr[dir];
                if(x+dx < 0 || x+dx >= w) continue;
                if(y+dy < 0 || y+dy >= h) continue;
                
                ztmp = imgrasters->zbuffer[(y+dy)*w + (x+dx)];
                if(ztmp==1) continue;
                ztmp  = ztmp * (sc->cam->maxdist-sc->cam->mindist) + sc->cam->mindist;
                rtmp = sc->cam->castRay(x, y, w, h);
                vtmp = rtmp.org + rtmp.dir*ztmp;
                dv = vtmp-v;
                d = dv.len();
                ao += fabs(dot(n, dv))*(1.0/(1.0+d))/d; //divide by d at the end to normalize dot product
                nsamples++;
            }
            ao /= fmax(1, nsamples);
            renderimg->setPixel(x, y, colorToRGB(RGBToColor(renderimg->pixel(x, y)) * color(1,1,1)*clamp(ao*2, 0, 1)));
        }
    }
    
}

void RenderArea::rayTraceUnthreaded(){
    num_rays_traced = 0;
    int i, j, x, y, xmax, ymax, w=width(), h=height();
    int tilenum=0, totaltiles = ceil((double)w/tilesize) * ceil((double)h/tilesize);
    int colrgb;
    vertex col;
    ray r;
    
    clock_t begin = clock();
    for(i=0; i<w; i+=tilesize){
        xmax = std::min(w, i+tilesize);
        for(j=0; j<h; j+=tilesize){
            ymax = std::min(h, j+tilesize);
            //for each tile
            for(x=i; x<xmax; x++){
                for(y=j; y<ymax; y++){
                    r = sc->cam->castRay(x, y, w, h);
                    if(amboc){
                        double occamount = ambientOcclusion(r, sc->kdt);
                        occamount = clamp(occamount, 0, 0.5)*2;
                        col.set(occamount, occamount, occamount);
                    }
                    else col = traceRay(r, 1, sc);
                    colrgb = colorToRGB(col);
                    renderimg->setPixel(x, y, colrgb);
                }
            }
            tilenum++;
            if(tilenum%100==0) std::cout << "tile " << tilenum << " of " << totaltiles  << " (" << (int)((100.0*tilenum)/totaltiles) << "%)" << std::endl;
        }
    }
    std::cout << "all tiles done. " << std::endl;
    clock_t end = clock();
    double time_elapsed = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << num_rays_traced << " rays traced in " << time_elapsed << " seconds, or " << num_rays_traced/time_elapsed << " rays per second." << std::endl;
}

void RenderArea::pathTraceUnthreaded(){
    num_rays_traced = 0;
    clock_t begin = clock();
    
    if(pathTracingSamples == 0) return;
    int i, j, x, y, xmax, ymax, w=width(), h=height(), s;
    int tilenum=0, totaltiles = ceil((double)w/tilesize) * ceil((double)h/tilesize);
    color totalcol;
    int colrgb;
    ray r;
    
    for(i=0; i<w; i+=tilesize){
        xmax = std::min(w, i+tilesize);
        for(j=0; j<h; j+=tilesize){
            ymax = std::min(h, j+tilesize);
            
            //for each tile
            for(x=i; x<xmax; x++){
                for(y=j; y<ymax; y++){
                    for(s=1, totalcol=color(); s<=pathTracingSamples; s++){
                        r = sc->cam->castRay(x, y, w, h);
                        totalcol += tracePath(r, 1, sc);
                    }
                    colrgb = colorToRGB(totalcol*(1.0/pathTracingSamples));
                    renderimg->setPixel(x, y, colrgb);
                }
            }
            tilenum++;
            if(tilenum%100==0) std::cout << "tile " << tilenum << " of " << totaltiles  << " (" << (int)((100.0*tilenum)/totaltiles) << "%)" << std::endl;
        }
    }
    
    std::cout << "all tiles done. " << std::endl;
    clock_t end = clock();
    double time_elapsed = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << num_rays_traced << " rays traced in " << time_elapsed << " seconds, or " << num_rays_traced/time_elapsed << " rays per second." << std::endl;
}

void RenderArea::drawStereoGram(){
    int w = width(), h=height(), x, y;
    double z;
    int r, color, shift;
    sc->cam->maxdist = 10;
    generate_maps_vector(1, imgrasters, sc);
    int tilesize = 100;
    int *pattern = new int[tilesize*tilesize];
    for(int i=0; i<tilesize*tilesize; i++){
            pattern[i] = rand();
    }
    
    for(y=0; y<h; y++){
        for(x=0; x<w; x++){
            z = imgrasters->zbuffer[y*w + x];
            shift = 31.0 - z*31.0;
            color = pattern[(y%tilesize) * tilesize + (x%(tilesize-shift))];
            renderimg->setPixel(x, y, color);
        }
    }
    
    delete pattern;
}

QColor colorToQColor(const color& c){
    return QColor((int)clamp(c.r*255, 0, 255), (int)clamp(c.g*255, 0, 255), (int)clamp(c.b*255, 0, 255));
}