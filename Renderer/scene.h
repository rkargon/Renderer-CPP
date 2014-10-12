//
//  scene.h
//  Renderer
//
//  Created by Raphael Kargon on 6/8/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__scene__
#define __Renderer__scene__

#include <QtGui/QImage>
#include "camera.h"
#include "geom.h"
#include "kdtree.h"

enum tex_projection_t {TEX_PROJ_SPHERICAL, TEX_PROJ_CUBIC, TEXT_PROJ_CYLINDRICAL};

class material{
public:
    color diff_col;
    real diff_intensity;
    
    color spec_col;
    real spec_intensity;
    real spec_hardness;
    
    real refl_intensity;
    real alpha;
    real ior;
    
    QImage *col_tex;
    QImage *spec_tex;
    QImage *norm_tex;
    
    material(const color& dc=color(1,1,1), const real di=1, const color& sc=color(1,1,1), const real si=1, const real sh=128, const real ri=0, const real a=1, const real indxofrefr=1.3, QImage *c=nullptr, QImage *s=nullptr, QImage *n=nullptr);
    color getColor(real txu, real txv);
    vertex getNormal(real txu, real txv);
    color getSpecCol(real txu, real txv);
};

class lamp{
public:
    real intensity;
    int falloff;
    vertex loc;
    color col;
    
    lamp();
    lamp(const real i, const int falloff, const vertex& l, const color& c);
};

class world{
public:
    color horizoncol, zenithcol;
    bool isFlat, isSky;
    world(const color& hc = color(1,1,1), const color& zc=color(0.8,0.8,1), const bool flat=false, const bool sky=false);
    virtual color getColor(const ray& r);
};

//uses Rayleigh and Mie scattering to simulate the sky
class sky : public world{
public:
    static const int samples=8;
    static const int samples_lightray=6;
    
    vertex betaR, betaM; //Rayleigh and Mie scattering coefficients at sea level
    real Hr, Hm; //Rayleigh and Mie scale heights
    real radiusEarth, radiusAtmo; //Earth and atmosphere radii
    vertex sundirection;
    real sunintensity;
    real g; //mean cosine
    
    sky(vertex betaR = vertex(5.5e-6, 13.0e-6, 22.4e-6), vertex betaM = vertex(21e-6, 21e-6, 21e-6), real Hr = 7994, real Hm = 1200, real radiusEarth=6360e3, real radiusAtmo=6420e3, vertex sundirection=vertex(0,1,0.2), real sunintensity=20, real g = 0.76);
    color getColor(const ray& r);
};

class scene {
public:
    camera *cam;
    std::vector<lamp*> lamps;
    world *w;
    mesh *obj;
    kdtree *kdt;
    
    scene(camera *c, std::vector<lamp*> l, world* wor, mesh *obj);
    scene();
};

#endif /* defined(__Renderer__material__) */
