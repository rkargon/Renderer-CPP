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
    zbuffer = new real[h*w];
    normalmap = new int[h*w];
    //set up render image
    renderimg = new QImage(w, h, QImage::Format_RGB32);
    
    //status label
    statuslbl = new QLabel("Raph Renderer 2014");
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
    clock_t begin = clock();
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
    //clock_t end = clock();
    //std::cout << real(end-begin)/CLOCKS_PER_SEC << std::endl;
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
    real zoomfactor = pow(1.005, -event->delta());
    sc->cam->zoom((real)zoomfactor);
    updateImage();
    repaint();
}

void RenderArea::resizeEvent(QResizeEvent *event){
    delete zbuffer;
    delete normalmap;
    delete renderimg;
    
    int h = height(), w = width();
    zbuffer = new real[h*w];
    normalmap = new int[h*w];
    renderimg = new QImage(w,h,QImage::Format_RGB32);
    updateImage();
}

/* DRAWING FUNCTIONS */

void RenderArea::drawWireFrame(){
    renderimg->fill(0xffffff);
    QPainter painter(renderimg);
    int w = width(), h = height();
    point2D<real> p1,p2;
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

void RenderArea::zBufferDraw(){
    renderimg->fill(0xffffff);
    int w = width(), h = height();
    real z, z1, z2, z3, dz21, dz31;
    int minx, miny, maxx, maxy;
    int A12, A23, A31, B12, B23, B31;
    int w0, w1, w2, w3, w1_row, w2_row, w3_row;
    int wsgn;
    vertex fcenter;
    point2D<real> p1, p2, p3, p;
    point2D<int> p1int, p2int, p3int, pint;
    color col;
    uint colrgb;
    
    //set up z buffer
    for(int i=0; i<w*h; i++) zbuffer[i] = 1;
    
    for(mesh *obj: sc->objects){
        for(face *f : obj->faces){
            //get pixels of vertices
            p1 = sc->cam->projectVertex(*f->vertices[0], w, h);
            p2 = sc->cam->projectVertex(*f->vertices[1], w, h);
            p3 = sc->cam->projectVertex(*f->vertices[2], w, h);
            if(isnan(p1.x)||isnan(p1.y) || isnan(p2.x)||isnan(p2.y) || isnan(p3.x)||isnan(p3.y)) continue;
            p1int.x=(int)p1.x;
            p1int.y=(int)p1.y;
            p2int.x=(int)p2.x;
            p2int.y=(int)p2.y;
            p3int.x=(int)p3.x;
            p3int.y=(int)p3.y;
            
            //z values = (z-min)/(max-min)
            z1 = (sc->cam->vertexDepth(*f->vertices[0])-sc->cam->mindist)/(sc->cam->maxdist-sc->cam->mindist);
            z2 = (sc->cam->vertexDepth(*f->vertices[1])-sc->cam->mindist)/(sc->cam->maxdist-sc->cam->mindist);
            z3 = (sc->cam->vertexDepth(*f->vertices[2])-sc->cam->mindist)/(sc->cam->maxdist-sc->cam->mindist);
            
            //store difference values. Makes interpolation later on slightly faster
            dz21 = z2-z1;
            dz31 = z3-z1;
            
            fcenter = f->center();
            col = calcLighting(fcenter, f->normal, *f->obj->mat, sc);
            colrgb = colorToRGB(col);
            
            //TODO smooth shading
            
            //triangle bounding box
            minx = std::min(std::min(p1int.x, p2int.x), p3int.x);
            miny = std::min(std::min(p1int.y, p2int.y), p3int.y);
            maxx = std::max(std::max(p1int.x, p2int.x), p3int.x);
            maxy = std::max(std::max(p1int.y, p2int.y), p3int.y);
            if(minx > w-1 || maxx<0 || miny>h-1 || maxy<0) continue; //face is off screen
            
            //clipint to screen
            minx = std::max(minx, 0);
            maxx = std::min(maxx, w-1);
            miny = std::max(miny, 0);
            maxy = std::min(maxy, h-1);
            
            
            //triangle edge setupint
            A12 = p1int.y-p2int.y;
            A23 = p2int.y-p3int.y;
            A31 = p3int.y-p1int.y;
            B12 = p2int.x-p1int.x;
            B23 = p3int.x-p2int.x;
            B31 = p1int.x-p3int.x;
            
            //initial barycentric coordinates at corner
            pint = point2D<int>(minx, miny);
            w1_row = orient2D(p2int,p3int,pint);
            w2_row = orient2D(p3int,p1int,pint);
            w3_row = orient2D(p1int,p2int,pint);
            w0 = orient2D(p1int, p2int, p3int);
            if(w0==0) continue;
            wsgn = signum(w0);
            //rasterize
            for(pint.y=miny; pint.y<=maxy; pint.y++){
                w1=w1_row;
                w2=w2_row;
                w3=w3_row;
                for(pint.x=minx; pint.x<=maxx; pint.x++){
                    if((signum(w1) == wsgn || !w1) && (signum(w2) == wsgn || !w2) && (signum(w3) == wsgn || !w3)){
                        //interpolate z value
                        z = z1 + w2*dz21/w0 + w3*dz31/w0;
                        if(z < zbuffer[w*pint.y+pint.x]){
                            zbuffer[w*pint.y+pint.x] = z;
                            renderimg->setPixel(pint.x, pint.y, colrgb);
                        }
                    }
                    w1+=A23;
                    w2+=A31;
                    w3+=A12;
                }
                w1_row += B23;
                w2_row += B31;
                w3_row += B12;
            }
        }
    }
}

void RenderArea::zBufferDraw_vector(){
    generate_maps_vector(5);
}

void RenderArea::paintNormalMap(){
    generate_maps_vector(3);
    renderimg->fill(0xffffff);
    for(int y=0; y<height(); y++){
        for(int x=0; x<width(); x++){
            renderimg->setPixel(x, y, normalmap[y*width()+x]);
        }
    }
}

//generates color, normal, and depth maps
// mapflags bit flags:
// 1 - depth map - (will always be generated anyway, needed for other maps)
// 2 - normal map
// 4 - color map
void RenderArea::generate_maps_vector(int mapflags){
    int i;
    int w = width(), h = height();
    int minx, miny, maxx, maxy;
    __v4sf z, z1, z2, z3, dz21, dz31; //interpolated z values
    __v4si w0, w1, w2, w3, w1_row, w2_row, w3_row, wsgn;
    __v4si pxmask; //whether each pixel is inside the triangle
    vertex fcenter;
    point2D<real> p1, p2, p3;
    point2D<int> p1int, p2int, p3int, pint;
    color col, col2, col3;
    vertex norm, norm2, norm3, normtmp;
    uint colrgb;
    std::map<meshvertex*, color> vertexcols;
    
    for(i=0; i<w*h; i++) zbuffer[i] = 1;
    if(mapflags & 2) for(i=0; i<w*h; i++) normalmap[i] = 0xffffff;
    if(mapflags & 4) renderimg->fill(0xffffff);
    
    for(mesh *obj: sc->objects){
        for(face *f : obj->faces){
            //get pixels of vertices
            p1 = sc->cam->projectVertex(*f->vertices[0], w, h);
            p2 = sc->cam->projectVertex(*f->vertices[1], w, h);
            p3 = sc->cam->projectVertex(*f->vertices[2], w, h);
            
            if(isnan(p1.x)||isnan(p1.y) || isnan(p2.x)||isnan(p2.y) || isnan(p3.x)||isnan(p3.y)) continue;
            p1int.x=(int)p1.x;
            p1int.y=(int)p1.y;
            p2int.x=(int)p2.x;
            p2int.y=(int)p2.y;
            p3int.x=(int)p3.x;
            p3int.y=(int)p3.y;
            
            //z values
            z1 = _mm_set1_ps((sc->cam->vertexDepth(*f->vertices[0])-sc->cam->mindist)/(sc->cam->maxdist-sc->cam->mindist));
            z2 = _mm_set1_ps((sc->cam->vertexDepth(*f->vertices[1])-sc->cam->mindist)/(sc->cam->maxdist-sc->cam->mindist));
            z3 = _mm_set1_ps((sc->cam->vertexDepth(*f->vertices[2])-sc->cam->mindist)/(sc->cam->maxdist-sc->cam->mindist));
            
            //store difference values. Makes interpolation later on slightly faster
            dz21 = z2-z1;
            dz31 = z3-z1;
            
            if (mapflags & 4){
                fcenter = f->center();
                colrgb = 1<<24; //unitialized color, largest byte is non-zero
            }
            
            //smooth shading, get vertex colors
            if((mapflags & (4+2)) && f->obj->smooth){
                norm = f->vertices[0]->vertexNormal();
                norm2 = f->vertices[1]->vertexNormal();
                norm3 = f->vertices[2]->vertexNormal();
                col = calcLighting(*f->vertices[0], norm, *f->obj->mat, sc);
                col2 = calcLighting(*f->vertices[1], norm2, *f->obj->mat, sc);
                col3 = calcLighting(*f->vertices[2], norm3, *f->obj->mat, sc);
            }
            
            //triangle bounding box
            minx = std::min(std::min(p1int.x, p2int.x), p3int.x);
            miny = std::min(std::min(p1int.y, p2int.y), p3int.y);
            maxx = std::max(std::max(p1int.x, p2int.x), p3int.x);
            maxy = std::max(std::max(p1int.y, p2int.y), p3int.y);
            if(minx > w-1 || maxx<0 || miny>h-1 || maxy<0) continue; //face is off screen
            
            //clip to screen
            minx = std::max(minx, 0);
            maxx = std::min(maxx, w-1);
            miny = std::max(miny, 0);
            maxy = std::min(maxy, h-1);
            
            
            //triangle edge setup
            pint = point2D<int>(minx, miny);
            EdgeVect e12, e23, e31;
            
            //initial barycentric coordinates at corner
            w1_row = e23.init(p2int, p3int, pint);
            w2_row = e31.init(p3int, p1int, pint);
            w3_row = e12.init(p1int, p2int, pint);
            w0 = _mm_set1_epi32(orient2D(p1int, p2int, p3int));
            //if w0 is zero continue
            if(_mm_movemask_epi8(_mm_cmpeq_epi32(w0,zeroveci))>0) { continue; }
            wsgn = _mm_cmpgt_epi32(w0, zeroveci);
            
            //rasterize
            for(pint.y=miny; pint.y<=maxy; pint.y+=EdgeVect::stepY){
                w1=w1_row;
                w2=w2_row;
                w3=w3_row;
                for(pint.x=minx; pint.x<=maxx; pint.x+=EdgeVect::stepX){
                    //each item in pxmask vector is >0 iff corresponding pixel should be drawn
                    //check if w1,2,3 have the same sign as w0, or are 0 themselves
                    pxmask = _mm_or_si128(_mm_cmpeq_epi32(_mm_cmpgt_epi32(w1, zeroveci), wsgn), _mm_cmpeq_epi32(w1,zeroveci));
                    pxmask = _mm_and_si128(pxmask, _mm_or_si128(_mm_cmpeq_epi32(_mm_cmpgt_epi32(w2, zeroveci), wsgn), _mm_cmpeq_epi32(w2,zeroveci)));
                    pxmask = _mm_and_si128(pxmask, _mm_or_si128(_mm_cmpeq_epi32(_mm_cmpgt_epi32(w3, zeroveci), wsgn), _mm_cmpeq_epi32(w3,zeroveci)));
                    //interpolate z value
                    //z = z1 + w2*dz21/w0 + w3*dz31/w0
                    z = z1 + _mm_div_ps(_mm_mul_ps(_mm_cvtepi32_ps(w2),dz21),_mm_cvtepi32_ps(w0)) + _mm_div_ps(_mm_mul_ps(_mm_cvtepi32_ps(w3),dz31),_mm_cvtepi32_ps(w0));
                    
                    for(i=0; i<4; i++){
                        if(pint.x+i<w && pxmask[i]!=0 && z[i] < zbuffer[w*pint.y+pint.x+i]){
                            if(mapflags & 2){
                                if(f->obj->smooth) normtmp = lerp(norm, norm2, norm3, real(w1[i])/w0[i], real(w2[i])/w0[i], real(w3[i])/w0[i]);
                                else normtmp = f->normal;
                                normalmap[pint.y*w + pint.x+i] = normalToRGB(normtmp);
                            }
                            if(mapflags & 4){
                                if(f->obj->smooth) colrgb = colorToRGB(lerp(col, col2, col3, real(w1[i])/w0[i], real(w2[i])/w0[i], real(w3[i])/w0[i]));
                                else if(colrgb>>24) colrgb = colorToRGB(calcLighting(fcenter, f->normal, *f->obj->mat, sc));
                                renderimg->setPixel(pint.x+i, pint.y, colrgb);
                            }
                            zbuffer[w*pint.y+pint.x+i] = z[i];
                        }
                    }
                    
                    w1+=e23.oneStepX;
                    w2+=e31.oneStepX;
                    w3+=e12.oneStepX;
                }
                w1_row += e23.oneStepY;
                w2_row += e31.oneStepY;
                w3_row += e12.oneStepY;
            }
        }
    }
}

//not really SSAO, should probably fix some tweaks
void RenderArea::SSAO()
{
    int w=width(), h=height();
    generate_maps_vector(7); //generate depth and normal maps
    //renderimg->fill(0xffffff);
    real ao, d, z, ztmp;
    vertex v, dv, vtmp, n;
    ray r, rtmp;
    int x, y, dy, dx, dir, nsamples;
    int dx_arr[4] = {-3, 3, 0, 0};
    int dy_arr[4] = {0, 0, -3, 3};
    
    for(y=0; y<h; y++){
        for(x=0; x<w; x++){
            z = zbuffer[y*w + x];
            if(z==1) continue;
            z = z * (sc->cam->maxdist-sc->cam->mindist) + sc->cam->mindist;
            r = sc->cam->castRay(x, y, w, h);
            v = r.org + r.dir*z;
            n = RGBToNormal(normalmap[y*w+x]);
            
            ao = 0;
            nsamples=0;
            for(dir=0; dir<4; dir++){
                dx = dx_arr[dir];
                dy = dy_arr[dir];
                if(x+dx < 0 || x+dx >= w) continue;
                if(y+dy < 0 || y+dy >= h) continue;
                
                ztmp = zbuffer[(y+dy)*w + (x+dx)];
                if(ztmp==1) continue;
                ztmp  = ztmp * (sc->cam->maxdist-sc->cam->mindist) + sc->cam->mindist;
                rtmp = sc->cam->castRay(x, y, w, h);
                vtmp = rtmp.org + rtmp.dir*ztmp;
                dv = vtmp-v;
                d = dv.len();
                ao += _abs(dot(n, dv))*(1.0/(1.0+d))/d; //divide by d at the end to normalize dot product
                nsamples++;
            }
            ao /= _max(1, nsamples);
            renderimg->setPixel(x, y, colorToRGB(RGBToColor(renderimg->pixel(x, y)) * color(1,1,1)*clamp(ao*2, 0, 1)));
        }
    }
    
}

void RenderArea::rayTraceUnthreaded(){
    num_rays_traced = 0;
    int i, j, x, y, xmax, ymax, w=width(), h=height();
    int tilenum=0, totaltiles = ceil((real)w/tilesize) * ceil((real)h/tilesize);
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
                        real occamount = ambientOcclusion(r, sc->kdt);
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
    real time_elapsed = real(end - begin) / CLOCKS_PER_SEC;
    std::cout << num_rays_traced << " rays traced in " << time_elapsed << " seconds, or " << num_rays_traced/time_elapsed << " rays per second." << std::endl;
}

void RenderArea::pathTraceUnthreaded(){
    num_rays_traced = 0;
    clock_t begin = clock();
    
    if(pathTracingSamples == 0) return;
    int i, j, x, y, xmax, ymax, w=width(), h=height(), s;
    int tilenum=0, totaltiles = ceil((real)w/tilesize) * ceil((real)h/tilesize);
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
    real time_elapsed = real(end - begin) / CLOCKS_PER_SEC;
    std::cout << num_rays_traced << " rays traced in " << time_elapsed << " seconds, or " << num_rays_traced/time_elapsed << " rays per second." << std::endl;
}

void RenderArea::drawStereoGram(){
    int w = width(), h=height(), x, y;
    real z;
    int r, color, shift;
    sc->cam->maxdist = 10;
    generate_maps_vector(1);
    int tilesize = 100;
    int *pattern = new int[tilesize*tilesize];
    for(int i=0; i<tilesize*tilesize; i++){
            pattern[i] = rand();
    }
    
    for(y=0; y<h; y++){
        for(x=0; x<w; x++){
            z = zbuffer[y*w + x];
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

__v4si EdgeVect::init(const point2D<int> &v0, const point2D<int> &v1, const point2D<int> &origin){
    // Edge setup
    int A = v0.y - v1.y, B = v1.x - v0.x;
    int C = v0.x*v1.y - v0.y*v1.x;
    
    //step deltas
    oneStepX = _mm_set1_epi32(A*stepX);
    oneStepY = _mm_set1_epi32(B*stepY);
    
    __v4si x = _mm_set1_epi32(origin.x) + _mm_set_epi32(3,2,1,0);
    __v4si y = _mm_set1_epi32(origin.y);
    
    //barycentric coordinates at edges:
    __v4si out = muli32(_mm_set1_epi32(A),x) + muli32(_mm_set1_epi32(B),y) + _mm_set1_epi32(C);
    return out;
}