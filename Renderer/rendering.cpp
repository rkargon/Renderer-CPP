//
//  rendering.cpp - the magic happens here!
//  Renderer
//
//  Created by Raphael Kargon on 6/8/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "rendering.h"

int num_rays_traced = 0;

//PRECONDITION: N is normalized
color calcLighting(const vertex& v, const vertex& n, const material& mat, scene* sc){
    color lightcol, speclightcol;
    for(lamp *l : sc->lamps){
        vertex lampvect = l->loc - v;
        vertex lampvnorm = lampvect.unitvect();
        real dotprod = dot(n,lampvnorm);
        vertex view = sc->cam->viewVector(v).unitvect();
        //if lamp and view are on different sides of face, then one is looking at underside of face)
        if(dotprod * dot(n, view) > 0) continue;
        real dstsqr = lampvect.lensqr();
        if(dstsqr==0){
            if(l->col.r) lightcol.r++;
            if(l->col.g) lightcol.g++;
            if(l->col.b) lightcol.b++;
            continue;
        }
        //Phong specular reflection
        vertex refl = lampvnorm.reflection(n);
        real spec_intensity = l->intensity * mat.spec_intensity * _max(0, _pow(dot(view, refl), mat.spec_hardness));
        real diff_intensity = l->intensity * _abs(dotprod)*mat.diff_intensity;
        //calc falloff
        diff_intensity /= dstsqr;
        spec_intensity /= dstsqr;
        lightcol += l->col*diff_intensity;
        speclightcol += l->col*spec_intensity;
    }
    color col = mat.diff_col * lightcol;
    color spcol = mat.spec_col * speclightcol;
    return col+spcol;
}

color traceRay(const ray& viewray, int depth, scene* sc){
    num_rays_traced++;
    vertex tuv;
    face *f = kdtree::rayTreeIntersect(sc->kdt, viewray, false, &tuv);
    if(f==nullptr) return sc->w->getColor(viewray);
    else{
        vertex v = viewray.org + viewray.dir*tuv.t; //calculate vertex location
        
        //interpolate texture coordinates
        real txu,txv;
        txu = f->vertices[0]->tex_u*(1-tuv.u-tuv.v) + f->vertices[1]->tex_u*tuv.u + f->vertices[2]->tex_u*tuv.v;
        txv = f->vertices[0]->tex_v*(1-tuv.u-tuv.v) + f->vertices[1]->tex_v*tuv.u + f->vertices[2]->tex_v*tuv.v;
        
        //calculate normal
        vertex n;
        if(f->obj->smooth){
            n = (f->vertices[0]->vertexNormal())*(1-tuv.u-tuv.v) + (f->vertices[1]->vertexNormal())*tuv.u + (f->vertices[2]->vertexNormal())*tuv.v;
            n.normalize();
        }
        else n = f->normal;
        //n = sc->obj->mat->getNormal(txu, txv);
        
        real ndotray = dot(n, viewray.dir);
        
        color lightcol, speclightcol;
        for(lamp *l : sc->lamps){
            vertex lampvect = l->loc - v;
            vertex lampvnorm = lampvect.unitvect();
            if(dot(n, lampvect)*ndotray>0) continue;
            ray lampray(v, lampvnorm);
            if(kdtree::rayTreeIntersect(sc->kdt, lampray, true, &tuv)!=nullptr) continue;
            real dotprod = dot(lampvnorm, n);
            real dstsqr = lampvect.lensqr();
            if(dstsqr==0){
                //if lamp is on the face's center, add full brightness to each color that the lamp emits
                if(l->col.r) lightcol.r++;
                if(l->col.g) lightcol.g++;
                if(l->col.b) lightcol.b++;
                continue;
            }
            //Phong shading
            vertex lamprefl = lampvnorm.reflection(n);
            real spec_intensity = l->intensity * f->obj->mat->spec_intensity * _max(0, _pow(dot(viewray.dir, lamprefl), f->obj->mat->spec_hardness));
            real diff_intensity = l->intensity *_abs(dotprod) * f->obj->mat->diff_intensity;
            //calculate falloff
            switch(l->falloff){
                case 2:
                    diff_intensity /= dstsqr;
                    spec_intensity /= dstsqr;
                    break;
                case 1:
					diff_intensity /= _sqrt(dstsqr);
					spec_intensity /= _sqrt(dstsqr);
                    break;
                case 0:
                    break;
            }
            lightcol += l->col * diff_intensity;
            speclightcol += l->col * spec_intensity;
        }
        
        color totcol = f->obj->mat->getColor(txu, txv) * lightcol + f->obj->mat->getSpecCol(txu, txv) * speclightcol;
        
        //an approximation. Assumes normals point outside and doesn't really deal with with concentric/intersecting objects.Currently assumes 'outside' of every object is air.
        //also doesn't do fresnel formula, reflection and refraction are handled separately, except for total internal reflection
        if(depth < RAY_DEPTH){
            vertex refl = viewray.dir.reflection(n);
            if(f->obj->mat->alpha < 1){
                real n1, n2;
                vertex transray;
                if(ndotray < 0){
                    n1 = 1;
                    n2 = f->obj->mat->ior;
                }
                else{
                    n1 = f->obj->mat->ior;
                    n2 = 1;
                }
                vertex raynorm = n * ndotray;
                vertex raytang = viewray.dir - raynorm;
                vertex transtang = raytang * (n1/n2);
                real transsinsquared = transtang.lensqr();
                if(transsinsquared > 1) transray = refl; //total reflection
                else{
                    vertex transnorm = n * signum(ndotray) *_sqrt(1-transsinsquared);
                    transray = transnorm + transtang;
                }
                color transcol = traceRay(ray(v, transray), depth+1, sc);
                totcol = lerp(transcol, totcol, f->obj->mat->alpha);
            }
            if(f->obj->mat->refl_intensity > 0){
                color refcol = traceRay(ray(v, refl), depth+1, sc);
                totcol = lerp(totcol, refcol, f->obj->mat->refl_intensity);
            }
        }
        return totcol;
    }
}

real ambientOcclusion(const ray& viewray, kdtree *kdt, int samples){
    //srand(0);
    vertex tuv;
    face *f = kdtree::rayTreeIntersect(kdt, viewray, false, &tuv);
    if(f==nullptr) return 1;
    else{
        vertex v = viewray.org + viewray.dir*tuv.t; //calculate vertex location
        //calculate normal
        vertex n;
        if(f->obj->smooth){
            n = (f->vertices[0]->vertexNormal())*(1-tuv.u-tuv.v) + (f->vertices[1]->vertexNormal())*tuv.u + (f->vertices[2]->vertexNormal())*tuv.v;
            n.normalize();
        }
        else n = f->normal;
        //real ndotray = dot(n, viewray.dir);
        
        real occ_amount = 0; //amount of ambient occlusion, ie how much current point is illuminated by the background.
        ray testray(v, vertex(0,0,0));
        for(int i=1; i<=samples; i++){
            testray.dir = randomDirection();
            if(dot(testray.dir, n) < 0) continue;
            else if (kdtree::rayTreeIntersect(kdt, testray, true, nullptr)==nullptr) occ_amount++;
        }
        occ_amount /= samples;
        return occ_amount;
    }
    srand(time(0));
}