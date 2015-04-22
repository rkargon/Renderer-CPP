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
#include "mesh.h"
#include "kdtree.h"

class material{
public:
    color diff_col; //The diffuse color of the material
    double diff_intensity; //The intensity of the diffuse color.
    
    color spec_col; //The color of specular highlights
    double spec_intensity; //The intensity of specular highlights
    double spec_hardness; //The hardness of specular highlights (according to Phong reflection model)
    
    double refl_intensity; //Amount of light reflected (0 = no relfection, 1 = mirror)
    double alpha; //Transparency of material. 0 = transparent, 1 = opaque
    double ior; //Index of refraction for transparency.
    
    //Textures for color, specular color, and normal map
    QImage *col_tex;
    QImage *spec_tex;
    QImage *norm_tex;
    
    material(const color& dc=color(1,1,1),
             const double di=1,
             const color& sc=color(1,1,1),
             const double si=1,
             const double sh=128,
             const double ri=0.3,
             const double a=1,
             const double indxofrefr=1.3,
             QImage *c=nullptr,
             QImage *s=nullptr,
             QImage *n=nullptr);
    
    //Get color, normal vector, and specular color using given texture coordinate
    color getColor(double txu, double txv);
    vertex getNormal(double txu, double txv);
    color getSpecCol(double txu, double txv);
};

class lamp{
public:
    double intensity;
    
    //Currently, falloff is not used and is assumed to be inverse quadratic
    int falloff; //The falloff exponent of the lamp's brightness
    //0 - constant falloff
    //1 - inverse linear falloff
    //2 - inverse quadratic falloff
    
    vertex loc;
    color col;
    
    lamp();
    lamp(const double i, const int falloff, const vertex& l, const color& c);
};

class world{
public:
    color horizoncol, zenithcol; //The color at horizon level, and the the color at the top of the sky
    bool isFlat; //If isFlat, all colors are horizoncol and zenithcol is ignored.
    world(const color& hc = color(1,1,1), const color& zc=color(0.8,0.8,1), const bool flat=false);
    virtual color getColor(const ray& r);
};

//uses Rayleigh and Mie scattering to simulate the sky
class sky : public world{
public:
    static const int samples=8;
    static const int samples_lightray=6;
    
    vertex betaR, betaM; //Rayleigh and Mie scattering coefficients at sea level
    double Hr, Hm; //Rayleigh and Mie scale heights
    double radiusEarth, radiusAtmo; //Earth and atmosphere radii
    vertex sundirection; //vector pointing in direction of the sun
    double sunintensity; //Intensity of the sun
    double g; //mean cosine
    
    sky(vertex betaR = vertex(5.5e-6, 13.0e-6, 22.4e-6),
        vertex betaM = vertex(21e-6, 21e-6, 21e-6),
        double Hr = 7994,
        double Hm = 1200,
        double radiusEarth=6360e3,
        double radiusAtmo=6420e3,
        vertex sundirection=vertex(0,-1,0.4).unitvect(),
        double sunintensity=20,
        double g = 0.76);
    color getColor(const ray& r);
};

class scene {
public:
    camera *cam;
    std::vector<lamp*> lamps;
    world *w;
    std::vector<mesh*> objects;
    kdtree *kdt;
    
    scene(camera *c, std::vector<lamp*> l, world* wor, std::vector<mesh*> objects);
    scene();
};

#endif /* defined(__Renderer__material__) */
