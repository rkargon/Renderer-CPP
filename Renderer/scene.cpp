//
//  scene.cpp
//  Renderer
//
//  Created by Raphael Kargon on 6/8/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "scene.h"

material::material(const color& dc, const double di, const color& sc, const double si, const double sh, const double ri, const double a, const double indxofrefr, QImage *c, QImage *s, QImage *n)
    :diff_col(dc),
    diff_intensity(di),
    spec_col(sc),
    spec_intensity(si),
    spec_hardness(sh),
    refl_intensity(ri),
    alpha(a),
    ior(indxofrefr),
    col_tex(c),
    spec_tex(s),
    norm_tex(n){}

//Either returns diff_col, or if col_tex is defined, the pixel corresponding to the passed texture coordinates
color material::getColor(double txu, double txv){
    if(col_tex==nullptr) return diff_col;
    else return RGBToColor(col_tex->pixel(txu*col_tex->width(), txv*col_tex->height()));
}
//TODO: Return normal relative to a passed, default normal
vertex material::getNormal(double txu, double txv){
    if(norm_tex==nullptr) return vertex();
    else return RGBToColor(norm_tex->pixel(txu*norm_tex->width(), txv*norm_tex->height()));
}
//Returns spec_col, or the corresponding pixel of spec_tex to determine specular color
color material::getSpecCol(double txu, double txv){
    if(spec_tex==nullptr) return spec_col;
    else return RGBToColor(spec_tex->pixel(txu*spec_tex->width(), txv*spec_tex->height()));
}

lamp::lamp()
    :intensity(1),
    falloff(2),
    loc(vertex()),
    col(color(1,1,1)){}

lamp::lamp(const double i, const int fall, const vertex& l, const color& c)
    :intensity(i),
    falloff(fall),
    loc(l),
    col(c){}


world::world(const color& hc, const color& zc, const bool flat)
    :horizoncol(hc),
    zenithcol(zc),
    isFlat(flat){}

//PRECONDITION: dir is normalized
//If isFlat, horizoncol is returned.
//Otherwise, z component of r.z is used to interpolate between horizoncol and zenithcol
color world::getColor(const ray &r){
    if(isFlat) return horizoncol;
    else{
        return lerp(horizoncol, zenithcol, fabs(r.dir.z)); //not really a linear progression from horizon to zenith, but close enough
    }
}

sky::sky(vertex betaR, vertex betaM, double Hr, double Hm, double radiusEarth, double radiusAtmo, vertex sundirection, double sunintensity, double g)
    :betaR(betaR),
    betaM(betaM),
    Hr(Hr),
    Hm(Hm),
    radiusEarth(radiusEarth),
    radiusAtmo(radiusAtmo),
    sundirection(sundirection.unitvect()),
    sunintensity(sunintensity),
    g(g){}


//Source: http://scratchapixel.com/lessons/3d-advanced-lessons/simulating-the-colors-of-the-sky/atmospheric-scattering/
//Given aa view ray, calculates the color of sky in that direction using atmospheric scattering.
color sky::getColor(const ray& r){
    ray r_new = r;
    r_new.org.z+=radiusEarth;
    //get atmosphere intersection point
    double t_atmo;
    if(!raySphereIntersect(r_new, radiusAtmo, t_atmo)) return color();
    double segmentLength = t_atmo / samples;
    vertex raysegment = r_new.dir*segmentLength, raysample = r_new.org+raysegment*0.5; //set up initial sample and delta for each sample along ray
    
    //get rayleigh and mie phase functions
    double mu = dot(r_new.dir, sundirection);
    double phaseR = 3*(1+mu*mu)/(16*M_PI);
    double phaseM = 3*(1-g*g)*(1+mu*mu)/(8*M_PI*(2+g*g)*pow(1+g*g-2*g*mu, 1.5));
    
    double opticalDepthR=0, opticalDepthM=0;
    color colR, colM; //resulting colors of each type of scattering
    
    //for each sample along ray
    for(int i=0; i<samples; i++, raysample+=raysegment){
        //get sample altitude
        double height = raysample.len() - radiusEarth;
        //update optical depth
        double hr = exp(-height/Hr)*segmentLength;
        double hm = exp(-height/Hm)*segmentLength;
        opticalDepthR += hr;
        opticalDepthM += hm;
        
        //calc amount of sunlight coming along ray at sample point
        ray lightray(raysample, sundirection);
        double lightray_t;
        //cast ray from sample to sun
        raySphereIntersect(lightray, radiusAtmo, lightray_t);
        double segmentLength_light = lightray_t/samples_lightray;
        vertex lightraysegment = lightray.dir * segmentLength_light,  lightraysample = lightray.org+lightraysegment*0.5;
        double opticalDepthLightR = 0, opticalDepthLightM = 0;

        //test several points along this ray
        int j;
        for(j=0; j<samples_lightray; j++, lightraysample+=lightraysegment){
            double heightlight = lightraysample.len() - radiusEarth;
            if(heightlight<0) break; //discard ray if it is in shadow of the earth
            //calc optical depth at each subsample
            opticalDepthLightR += exp(-heightlight/Hr)*segmentLength_light;
            opticalDepthLightM += exp(-heightlight/Hm)*segmentLength_light;
        }
        if(j==samples_lightray){ //ray is not in shadow of earth
            //calc light color based on attenuation
            vertex tau = betaR*(opticalDepthR+opticalDepthLightR) + betaM*1.1*(opticalDepthM+opticalDepthLightM);
            vertex attenuation(exp(-tau.x), exp(-tau.y), exp(-tau.z));
            colR += hr * attenuation;
            colM += hm * attenuation;
        }
        else return color();
    }
    return (colR*betaR*phaseR + colM*betaM*phaseM)*sunintensity;
}

scene::scene()
:cam(nullptr),
lamps(),
w(nullptr),
objects(std::vector<mesh*>()),
kdt(nullptr){}

// automatically builds kdtree for object
scene::scene(camera *c, std::vector<lamp*> l, world* wor, std::vector<mesh*> os)
:cam(c),
lamps(l),
w(wor),
objects(os){
    std::vector<face*> allfaces;
    for(mesh *o : objects){
        allfaces.insert(allfaces.end(), o->faces.begin(), o->faces.end());
    }
    this->kdt = kdtree::buildTree(allfaces);
}
