//
//  scene.cpp
//  Renderer
//
//  Created by Raphael Kargon on 6/8/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "scene.h"

material::material(const color& dc, const real di, const color& sc, const real si, const real sh, const real ri, const real a, const real indxofrefr, QImage *c, QImage *s, QImage *n)
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

color material::getColor(real txu, real txv){
    if(col_tex==nullptr) return diff_col;
    else return RGBToColor(col_tex->pixel(txu*col_tex->width(), txv*col_tex->height()));
}
vertex material::getNormal(real txu, real txv){
    if(norm_tex==nullptr) return vertex();
    else return RGBToColor(norm_tex->pixel(txu*norm_tex->width(), txv*norm_tex->height()));
}
color material::getSpecCol(real txu, real txv){
    if(spec_tex==nullptr) return spec_col;
    else return RGBToColor(spec_tex->pixel(txu*spec_tex->width(), txv*spec_tex->height()));
}

lamp::lamp()
    :intensity(1),
    falloff(2),
    loc(vertex()),
    col(color(1,1,1)){}

lamp::lamp(const real i, const int fall, const vertex& l, const color& c)
    :intensity(i),
    falloff(fall),
    loc(l),
    col(c){}


world::world(const color& hc, const color& zc, const bool flat, const bool sky)
    :horizoncol(hc),
    zenithcol(zc),
    isFlat(flat),
    isSky(sky){}

//PRECONDITION: dir is normalized
color world::getColor(const ray &r){
    if(isFlat) return horizoncol;
    else{
        return lerp(horizoncol, zenithcol, _abs(r.dir.z)); //not really a linear progression from horizon to zenith, but close enough
    }
}

sky::sky(vertex betaR, vertex betaM, real Hr, real Hm, real radiusEarth, real radiusAtmo, vertex sundirection, real sunintensity, real g)
    :betaR(betaR),
    betaM(betaM),
    Hr(Hr),
    Hm(Hm),
    radiusEarth(radiusEarth),
    radiusAtmo(radiusAtmo),
    sundirection(sundirection.unitvect()),
    sunintensity(sunintensity),
    g(g){}

color sky::getColor(const ray& r){
    ray r_new = r;
    r_new.org.z+=radiusEarth;
    //get atmosphere intersection point
    real t_atmo;
    if(!raySphereIntersect(r_new, radiusAtmo, t_atmo)) return color();
    real segmentLength = t_atmo / samples;
    vertex raysegment = r_new.dir*segmentLength, raysample = r_new.org+raysegment*0.5; //set up initial sample and delta for each sample along ray
    
    //get rayleigh and mie phase functions
    real mu = dot(r_new.dir, sundirection);
    real phaseR = 3*(1+mu*mu)/(16*M_PI);
    real phaseM = 3*(1-g*g)*(1+mu*mu)/(8*M_PI*(2+g*g)*_pow(1+g*g-2*g*mu, 1.5));
    
    real opticalDepthR=0, opticalDepthM=0;
    color colR, colM; //resulting colors of each type of scattering
    
    //for each sample along ray
    for(int i=0; i<samples; i++, raysample+=raysegment){
        //get sample altitude
        real height = raysample.len() - radiusEarth;
        //update optical depth
        real hr = _exp(-height/Hr)*segmentLength;
        real hm = _exp(-height/Hm)*segmentLength;
        opticalDepthR += hr;
        opticalDepthM += hm;
        
        //calc amount of sunlight coming along ray at sample point
        ray lightray(raysample, sundirection);
        real lightray_t;
        //cast ray from sample to sun
        raySphereIntersect(lightray, radiusAtmo, lightray_t);
        real segmentLength_light = lightray_t/samples_lightray;
        vertex lightraysegment = lightray.dir * segmentLength_light,  lightraysample = lightray.org+lightraysegment*0.5;
        real opticalDepthLightR = 0, opticalDepthLightM = 0;

        //test several points along this ray
        int j;
        for(j=0; j<samples_lightray; j++, lightraysample+=lightraysegment){
            real heightlight = lightraysample.len() - radiusEarth;
            if(heightlight<0) break; //discard ray if it is in shadow of the earth
            //calc optical depth at each subsample
            opticalDepthLightR += _exp(-heightlight/Hr)*segmentLength_light;
            opticalDepthLightM += _exp(-heightlight/Hm)*segmentLength_light;
        }
        if(j==samples_lightray){ //ray is not in shadow of earth
            //calc light color based on attenuation
            vertex tau = betaR*(opticalDepthR+opticalDepthLightR) + betaM*1.1*(opticalDepthM+opticalDepthLightM);
            vertex attenuation(_exp(-tau.x), _exp(-tau.y), _exp(-tau.z));
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
obj(nullptr),
kdt(nullptr){}

scene::scene(camera *c, std::vector<lamp*> l, world* wor, mesh *o)
:cam(c),
lamps(l),
w(wor),
obj(o),
kdt(nullptr){}
