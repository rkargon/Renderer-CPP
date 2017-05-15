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
color material::get_color(double txu, double txv){
    if(col_tex==nullptr) return diff_col;
    else return rgb_to_color(col_tex->pixel(txu*col_tex->width(), txv*col_tex->height()));
}
//TODO: Return normal relative to a passed, default normal
vertex material::get_normal(double txu, double txv){
    if(norm_tex==nullptr) return vertex();
    else return rgb_to_color(norm_tex->pixel(txu*norm_tex->width(), txv*norm_tex->height()));
}
//Returns spec_col, or the corresponding pixel of spec_tex to determine specular color
color material::get_spec_col(double txu, double txv){
    if(spec_tex==nullptr) return spec_col;
    else return rgb_to_color(spec_tex->pixel(txu*spec_tex->width(), txv*spec_tex->height()));
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
col(c),
is_sun(false){}

lamp::lamp(const double i, const int fall, const vertex& l, const color& c, const vertex& sun_direction)
:intensity(i),
falloff(fall),
loc(l),
col(c),
sun_direction(sun_direction),
is_sun(true){}

world::world(const color& hc, const color& zc, const bool flat, const double ambient_intensity)
    :horizon_col(hc),
    zenith_col(zc),
    is_flat(flat),
    ambient_intensity(ambient_intensity){}

//PRECONDITION: dir is normalized
//If is_flat, horizon_col is returned.
//Otherwise, z component of r.z is used to interpolate between horizon_col and zenith_col
color world::get_color(const ray &r){
    if(is_flat) return horizon_col;
    else{
        return lerp(horizon_col, zenith_col, fabs(r.dir.z)); //not really a linear progression from horizon to zenith, but close enough
    }
}

sky::sky(vertex beta_r, vertex beta_m, double Hr, double Hm, double radius_earth, double radius_atmo, vertex sun_direction, double sun_intensity, double g)
    :beta_r(beta_r),
    beta_m(beta_m),
    Hr(Hr),
    Hm(Hm),
    radius_earth(radius_earth),
    radius_atmo(radius_atmo),
    sun_direction(sun_direction.unitvect()),
    sun_intensity(sun_intensity),
    g(g){}


//Source: http://scratchapixel.com/lessons/3d-advanced-lessons/simulating-the-colors-of-the-sky/atmospheric-scattering/
//Given aa view ray, calculates the color of sky in that direction using atmospheric scattering.
color sky::get_color(const ray& r){
    ray r_new = r;
    r_new.org.z+=radius_earth;
    //get atmosphere intersection point
    double t_atmo;
    if(!ray_sphere_intersect(r_new, radius_atmo, t_atmo)) return color();
    double segment_length = t_atmo / samples;
    vertex ray_segment = r_new.dir*segment_length, ray_sample = r_new.org+ray_segment*0.5; //set up initial sample and delta for each sample along ray
    
    //get rayleigh and mie phase functions
    double mu = dot(r_new.dir, sun_direction);
    double phase_r = 3*(1+mu*mu)/(16*M_PI);
    double phase_m = 3*(1-g*g)*(1+mu*mu)/(8*M_PI*(2+g*g)*pow(1+g*g-2*g*mu, 1.5));
    
    double optical_depth_r=0, optical_depth_m=0;
    color col_r{}, col_m{}; //resulting colors of each type of scattering
    //for each sample along ray
    for(int i=0; i<samples; i++, ray_sample+=ray_segment){
        //get sample altitude
        double height = ray_sample.len() - radius_earth;
        //update optical depth
        double hr = exp(-height/Hr)*segment_length;
        double hm = exp(-height/Hm)*segment_length;
        optical_depth_r += hr;
        optical_depth_m += hm;
        
        //calc amount of sunlight coming along ray at sample point
        ray lightray(ray_sample, sun_direction);
        double lightray_t;
        //cast ray from sample to sun
        ray_sphere_intersect(lightray, radius_atmo, lightray_t);
        double segment_length_light = lightray_t/samples_lightray;
        vertex lightray_segment = lightray.dir * segment_length_light;
        vertex lightray_sample = lightray.org+lightray_segment*0.5;
        double optical_depth_light_r = 0, optical_depth_light_m = 0;

        //test several points along this ray
        int j;
        for(j=0; j<samples_lightray; j++, lightray_sample+=lightray_segment){
            double heightlight = lightray_sample.len() - radius_earth;
            if(heightlight<0) break; //discard ray if it is in shadow of the earth
            //calc optical depth at each subsample
            optical_depth_light_r += exp(-heightlight/Hr)*segment_length_light;
            optical_depth_light_m += exp(-heightlight/Hm)*segment_length_light;
        }
        if(j==samples_lightray){ //ray is not in shadow of earth
            //calc light color based on attenuation
            vertex tau = beta_r*(optical_depth_r+optical_depth_light_r) + beta_m*1.1*(optical_depth_m+optical_depth_light_m);
            vertex attenuation(exp(-tau.x), exp(-tau.y), exp(-tau.z));
            col_r += hr * attenuation;
            col_m += hm * attenuation;
        }
        else return color();
    }
    return (col_r*beta_r*phase_r + col_m*beta_m*phase_m)*sun_intensity;
}

scene::scene()
:cam(nullptr),
lamps(),
w(nullptr),
objects(std::vector<mesh*>()),
kdt(nullptr){}

// automatically builds kdtree for object
scene::scene(camera *c, std::vector<lamp*> l, world* wor, std::vector<mesh*> objects, distance_estimator *de_obj, material *de_mat)
:cam(c),
lamps(l),
w(wor),
objects(objects),
de_obj(de_obj),
de_mat(de_mat){
    std::vector<face*> allfaces;
    for(mesh *o : objects){
        allfaces.insert(allfaces.end(), o->faces.begin(), o->faces.end());
    }
    this->kdt = kdtree::build_tree(allfaces);
}
