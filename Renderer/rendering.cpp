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
// Calculate lighting for triangle rasterization.
// Given a vertex, a normal, and a material, will calculate the color of it, taking into account incident light from lamps.
// Does not take into account other geometry which could cause shadows, reflections etc.
color calcLighting(const vertex& v, const vertex& n, const material& mat, const scene* sc){
    color lightcol, speclightcol;
    for(lamp *l : sc->lamps){
        vertex lampvect = l->loc - v;
        vertex lampvnorm = lampvect.unitvect();
        double dotprod = dot(n,lampvnorm);
        vertex view = sc->cam->viewVector(v).unitvect();
        //if lamp and view are on different sides of face, then one is looking at underside of face)
        if(dotprod * dot(n, view) > 0) continue;
        double dstsqr = lampvect.lensqr();
        if(dstsqr==0){
            if(l->col.r) lightcol.r++;
            if(l->col.g) lightcol.g++;
            if(l->col.b) lightcol.b++;
            continue;
        }
        //Phong specular reflection
        vertex refl = lampvnorm.reflection(n);
        double spec_intensity = l->intensity * mat.spec_intensity * fmax(0, pow(dot(view, refl), mat.spec_hardness));
        double diff_intensity = l->intensity * fabs(dotprod)*mat.diff_intensity;
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

// Traces a ray from the camera through a scene, using ray-tracing, up to a certain depth
color traceRay(const ray& viewray, scene* sc, int depth){
    num_rays_traced++;
    vertex tuv;
    face *f = kdtree::rayTreeIntersect(sc->kdt, viewray, false, &tuv);
    if(f==nullptr) return sc->w->getColor(viewray);
    else{
        vertex v = viewray.org + viewray.dir*tuv.t; //calculate vertex location
        
        //interpolate texture coordinates
        double txu,txv;
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
        
        double ndotray = dot(n, viewray.dir);
        
        /* SPECULAR & DIFFUSE LIGHTING */
        color lightcol, speclightcol;
        for(lamp *l : sc->lamps){
            vertex lampvect = l->loc - v;
            vertex lampvnorm = lampvect.unitvect();
            if(dot(n, lampvect)*ndotray>0) continue; //make sure lamp is on same side of face as view
            ray lampray(v, lampvnorm);
            if(kdtree::rayTreeIntersect(sc->kdt, lampray, true, &tuv)!=nullptr) continue;
            double dotprod = dot(lampvnorm, n);
            double dstsqr = lampvect.lensqr();
            if(dstsqr==0){
                //if lamp is on the face's center, add full brightness to each color that the lamp emits
                if(l->col.r) lightcol.r++;
                if(l->col.g) lightcol.g++;
                if(l->col.b) lightcol.b++;
                continue;
            }
            //Phong shading
            vertex lamprefl = lampvnorm.reflection(n);
            double spec_intensity = l->intensity * f->obj->mat->spec_intensity * fmax(0, pow(dot(viewray.dir, lamprefl), f->obj->mat->spec_hardness));
            double diff_intensity = l->intensity *fabs(dotprod) * f->obj->mat->diff_intensity;
            //calculate falloff
            switch(l->falloff){
                case 2:
                    diff_intensity /= dstsqr;
                    spec_intensity /= dstsqr;
                    break;
                case 1:
					diff_intensity /= sqrt(dstsqr);
					spec_intensity /= sqrt(dstsqr);
                    break;
                case 0:
                    break;
            }
            lightcol += l->col * diff_intensity;
            speclightcol += l->col * spec_intensity;
        }
        
        color totcol = f->obj->mat->getColor(txu, txv) * lightcol + f->obj->mat->getSpecCol(txu, txv) * speclightcol;
        
        /* Ambient lighting from sky */
        if (sc->w->ambient_intensity > 0){
            ray normal_ray{v, n};
            if(kdtree::rayTreeIntersect(sc->kdt, normal_ray, true, nullptr)==nullptr){
                totcol += sc->w->ambient_intensity * sc->w->getColor(normal_ray);
            }
        }
        
        /* REFLECTION & REFRACTION */
        //an approximation. Assumes normals point outside and doesn't really deal with with concentric/intersecting objects.Currently assumes 'outside' of every object is air.
        //also doesn't do fresnel formula, reflection and refraction are handled separately, except for total internal reflection
        if(depth < RAY_DEPTH){
            vertex refl = viewray.dir.reflection(n);
            if(f->obj->mat->alpha < 1){
                double n1, n2;
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
                double transsinsquared = transtang.lensqr();
                if(transsinsquared > 1) transray = refl; //total reflection
                else{
                    vertex transnorm = n * signum(ndotray) *sqrt(1-transsinsquared);
                    transray = transnorm + transtang;
                }
                color transcol = traceRay(ray(v, transray), sc, depth+1);
                totcol = lerp(transcol, totcol, f->obj->mat->alpha);
            }
            if(f->obj->mat->refl_intensity > 0){
                color refcol = traceRay(ray(v, refl), sc, depth+1);
                totcol = lerp(totcol, refcol, f->obj->mat->refl_intensity);
            }
        }
        return totcol;
    }
}

// Traces a path through a scene, up to the given depth.
// Used for path-tracing, and calculates indirect lighting.
color tracePath(const ray& viewray, scene* sc, int depth){
    if(depth > RAY_DEPTH) return color();
    
    num_rays_traced++;
    vertex tuv;
    face *f = kdtree::rayTreeIntersect(sc->kdt, viewray, false, &tuv);
    if(f==nullptr) return sc->w->getColor(viewray);
    else{
        vertex v = viewray.org + viewray.dir*tuv.t; //calculate vertex location
    
        //calculate normal
        vertex n;
        if(f->obj->smooth){
            n = (f->vertices[0]->vertexNormal())*(1-tuv.u-tuv.v) + (f->vertices[1]->vertexNormal())*tuv.u + (f->vertices[2]->vertexNormal())*tuv.v;
            n.normalize();
        }
        else n = f->normal;
        
        //calculate incident light
        vertex inc_dir = f->obj->bsdf->getIncidentDirection(n, viewray.dir);
        color inc_col = tracePath(ray(v, inc_dir), sc, depth+1);
        //calculate returned light
        color return_col = f->obj->bsdf->getLight(inc_col, inc_dir, n, viewray.dir);
        return return_col;
    }
}

// Calculates ambient occlusion for a ray.
double ambientOcclusion(const ray& viewray, scene *sc){
    vertex tuv;
    face *f = kdtree::rayTreeIntersect(sc->kdt, viewray, false, &tuv);
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
        //double ndotray = dot(n, viewray.dir);
        
        double occ_amount = 0; //amount of ambient occlusion, ie how much current point is illuminated by the background.
        ray testray(v, vertex(0,0,0));
        for(int i=1; i<=AMB_OCC_SAMPLES; i++){
            testray.dir = randomDirection();
            if(dot(testray.dir, n) < 0) continue;
            else if (kdtree::rayTreeIntersect(sc->kdt, testray, true, nullptr)==nullptr) occ_amount++;
        }
        occ_amount /= AMB_OCC_SAMPLES;
        return occ_amount;
    }
}

color rayTraceDistanceField(const ray& viewray, scene *sc, int num_steps, int depth){
    num_rays_traced++;
    double t;
    int steps;
    if (!ray_march(viewray, *sc->de_obj, &t, &steps, num_steps)){
        return sc->w->getColor(viewray);
    } else {
        vertex v = viewray.org + t*viewray.dir;
        double ambient_occlusion = 1.0 - (steps / (double)num_steps);
        vertex n =  estimate_normal(v, *sc->de_obj);
        double ndotray = dot(n, viewray.dir);
        
        /* SPECULAR & DIFFUSE LIGHTING */
        color lightcol, speclightcol;
        for(lamp *l : sc->lamps){
            vertex lampvect = l->loc - v;
            vertex lampvnorm = lampvect.unitvect();
            if(dot(n, lampvect)*ndotray>0) continue; //make sure lamp is on same side of face as view
            ray lampray(v, lampvnorm);
            if(ray_march(lampray, *sc->de_obj, nullptr, nullptr, num_steps)){
                continue;
            }
            double dotprod = dot(lampvnorm, n);
            double dstsqr = lampvect.lensqr();
            if(dstsqr==0){
                //if lamp is on the face's center, add full brightness to each color that the lamp emits
                if(l->col.r) lightcol.r++;
                if(l->col.g) lightcol.g++;
                if(l->col.b) lightcol.b++;
                continue;
            }
            //Phong shading
            vertex lamprefl = lampvnorm.reflection(n);
            double spec_intensity = l->intensity * sc->de_mat->spec_intensity * fmax(0, pow(dot(viewray.dir, lamprefl), sc->de_mat->spec_hardness));
            double diff_intensity = l->intensity *fabs(dotprod) * sc->de_mat->diff_intensity;
            //calculate falloff
            switch(l->falloff){
                case 2:
                    diff_intensity /= dstsqr;
                    spec_intensity /= dstsqr;
                    break;
                case 1:
                    diff_intensity /= sqrt(dstsqr);
                    spec_intensity /= sqrt(dstsqr);
                    break;
                case 0:
                    break;
            }
            lightcol += l->col * diff_intensity;
            speclightcol += l->col * spec_intensity;
        }
        color totcol = sc->de_mat->diff_col * lightcol + sc->de_mat->spec_col * speclightcol;
        
        /* Ambient lighting from sky */
        if (sc->w->ambient_intensity > 0){
            ray normal_ray{v, n};
            // double t;
            if(!ray_march(normal_ray, *sc->de_obj, nullptr, nullptr, num_steps)){
                totcol += sc->w->ambient_intensity * sc->w->getColor(normal_ray) * ambient_occlusion;
            }
        }
        
        /* REFLECTION & REFRACTION */
        // An approximation. Assumes normals point outside and doesn't really deal with with concentric/intersecting objects.
        // Currently assumes 'outside' of every object is air.
        //also doesn't do fresnel formula, reflection and refraction are handled separately, except for total internal reflection
        if(depth < RAY_DEPTH){
            vertex refl = viewray.dir.reflection(n);
            if(sc->de_mat->alpha < 1){
                double n1, n2;
                vertex transray;
                if(ndotray < 0){
                    n1 = 1;
                    n2 = sc->de_mat->ior;
                }
                else{
                    n1 = sc->de_mat->ior;
                    n2 = 1;
                }
                vertex raynorm = n * ndotray;
                vertex raytang = viewray.dir - raynorm;
                vertex transtang = raytang * (n1/n2);
                double transsinsquared = transtang.lensqr();
                if(transsinsquared > 1) transray = refl; //total reflection
                else{
                    vertex transnorm = n * signum(ndotray) *sqrt(1-transsinsquared);
                    transray = transnorm + transtang;
                }
                // TODO
                color transcol = traceRay(ray(v, transray), sc, depth+1);
                totcol = lerp(transcol, totcol, sc->de_mat->alpha);
            }
            if(sc->de_mat->refl_intensity > 0){
                color refcol = traceRay(ray(v, refl), sc, depth+1);
                totcol = lerp(totcol, refcol, sc->de_mat->refl_intensity);
            }
        }
        return totcol;
    }
    
    vertex tuv;
    face *f = kdtree::rayTreeIntersect(sc->kdt, viewray, false, &tuv);
    if(f==nullptr) return sc->w->getColor(viewray);
    else{
        vertex v = viewray.org + viewray.dir*tuv.t; //calculate vertex location
        
        //interpolate texture coordinates
        double txu,txv;
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
        
        double ndotray = dot(n, viewray.dir);
        
        /* SPECULAR & DIFFUSE LIGHTING */
        color lightcol, speclightcol;
        for(lamp *l : sc->lamps){
            vertex lampvect = l->loc - v;
            vertex lampvnorm = lampvect.unitvect();
            if(dot(n, lampvect)*ndotray>0) continue; //make sure lamp is on same side of face as view
            ray lampray(v, lampvnorm);
            if(kdtree::rayTreeIntersect(sc->kdt, lampray, true, &tuv)!=nullptr) continue;
            double dotprod = dot(lampvnorm, n);
            double dstsqr = lampvect.lensqr();
            if(dstsqr==0){
                //if lamp is on the face's center, add full brightness to each color that the lamp emits
                if(l->col.r) lightcol.r++;
                if(l->col.g) lightcol.g++;
                if(l->col.b) lightcol.b++;
                continue;
            }
            //Phong shading
            vertex lamprefl = lampvnorm.reflection(n);
            double spec_intensity = l->intensity * f->obj->mat->spec_intensity * fmax(0, pow(dot(viewray.dir, lamprefl), f->obj->mat->spec_hardness));
            double diff_intensity = l->intensity *fabs(dotprod) * f->obj->mat->diff_intensity;
            //calculate falloff
            switch(l->falloff){
                case 2:
                    diff_intensity /= dstsqr;
                    spec_intensity /= dstsqr;
                    break;
                case 1:
                    diff_intensity /= sqrt(dstsqr);
                    spec_intensity /= sqrt(dstsqr);
                    break;
                case 0:
                    break;
            }
            lightcol += l->col * diff_intensity;
            speclightcol += l->col * spec_intensity;
        }
        
        color totcol = f->obj->mat->getColor(txu, txv) * lightcol + f->obj->mat->getSpecCol(txu, txv) * speclightcol;
        
        /* REFLECTION & REFRACTION */
        //an approximation. Assumes normals point outside and doesn't really deal with with concentric/intersecting objects.Currently assumes 'outside' of every object is air.
        //also doesn't do fresnel formula, reflection and refraction are handled separately, except for total internal reflection
        if(depth < RAY_DEPTH){
            vertex refl = viewray.dir.reflection(n);
            if(f->obj->mat->alpha < 1){
                double n1, n2;
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
                double transsinsquared = transtang.lensqr();
                if(transsinsquared > 1) transray = refl; //total reflection
                else{
                    vertex transnorm = n * signum(ndotray) *sqrt(1-transsinsquared);
                    transray = transnorm + transtang;
                }
                color transcol = traceRay(ray(v, transray), sc, depth+1);
                totcol = lerp(transcol, totcol, f->obj->mat->alpha);
            }
            if(f->obj->mat->refl_intensity > 0){
                color refcol = traceRay(ray(v, refl), sc, depth+1);
                totcol = lerp(totcol, refcol, f->obj->mat->refl_intensity);
            }
        }
        return totcol;
    }
}

/* Rasterization */

//generates color, normal, and depth maps
// mapflags bit flags:
// 1 - depth map - (will always be generated anyway, needed for other maps)
// 2 - normal map
// 4 - color map
//reference implementation, no vector operations
void generate_maps(int mapflags, raster *imgrasters, scene *sc){
    int w = imgrasters->width(), h = imgrasters->height();
    double z, z1, z2, z3, dz21, dz31;
    int minx, miny, maxx, maxy;
    int A12, A23, A31, B12, B23, B31;
    int w0, w1, w2, w3, w1_row, w2_row, w3_row;
    int wsgn;
    vertex fcenter;
    point2D<double> p1, p2, p3, p;
    point2D<int> p1int, p2int, p3int, pint;
    color col, col2, col3;
    vertex norm, norm2, norm3, normtmp;
    uint colrgb = 1;
    
    std::fill_n(imgrasters->zbuffer, w*h, 1);
    if(mapflags & 2) std::fill_n(imgrasters->normbuffer, w*h, 0xffffff);
    if(mapflags & 4) std::fill_n(imgrasters->colbuffer, w*h, 0xffffff);
    
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
            
            //clipint to screen
            minx = std::max(minx, 0);
            maxx = std::min(maxx, w-1);
            miny = std::max(miny, 0);
            maxy = std::min(maxy, h-1);
            
            
            //triangle edge setup
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
                        if(z < imgrasters->zbuffer[w*pint.y+pint.x]){
                            if(mapflags & 2){
                                if(f->obj->smooth) normtmp = lerp(norm, norm2, norm3, double(w1)/w0, double(w2)/w0, double(w3)/w0);
                                else normtmp = f->normal;
                                imgrasters->normbuffer[pint.y*w + pint.x] = normalToRGB(normtmp);
                            }
                            if(mapflags & 4){
                                if(f->obj->smooth) colrgb = colorToRGB(lerp(col, col2, col3, double(w1)/w0, double(w2)/w0, double(w3)/w0));
                                else if(colrgb>>24) colrgb = colorToRGB(calcLighting(fcenter, f->normal, *f->obj->mat, sc));
                                imgrasters->colbuffer[pint.y*w + pint.x] = colrgb;
                            }
                            imgrasters->zbuffer[w*pint.y + pint.x] = z;
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

//generates color, normal, and depth maps
// mapflags bit flags:
// 1 - depth map - (will always be generated anyway, needed for other maps)
// 2 - normal map
// 4 - color map
void generate_maps_vector(int mapflags, raster *imgrasters, scene *sc){
    int i;
    int w = imgrasters->width(), h = imgrasters->height();
    int minx, miny, maxx, maxy;
    __v4sf z, z1, z2, z3, dz21, dz31; //interpolated z values
    __v4si w0, w1, w2, w3, w1_row, w2_row, w3_row, wsgn;
    __v4si pxmask; //whether each pixel is inside the triangle
    vertex fcenter;
    point2D<double> p1, p2, p3;
    point2D<int> p1int, p2int, p3int, pint;
    color col, col2, col3;
    vertex norm, norm2, norm3, normtmp;
    uint colrgb = 1;
    
    std::fill_n(imgrasters->zbuffer, w*h, 1);
    if(mapflags & 2) std::fill_n(imgrasters->normbuffer, w*h, 0xffffff);
    if(mapflags & 4) std::fill_n(imgrasters->colbuffer, w*h, 0xffffff);
    
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
                        if(pint.x+i<w && pxmask[i]!=0 && z[i] < imgrasters->zbuffer[w*pint.y+pint.x+i]){
                            if(mapflags & 2){
                                if(f->obj->smooth) normtmp = lerp(norm, norm2, norm3, double(w1[i])/w0[i], double(w2[i])/w0[i], double(w3[i])/w0[i]);
                                else normtmp = f->normal;
                                imgrasters->normbuffer[pint.y*w + pint.x+i] = normalToRGB(normtmp);
                            }
                            if(mapflags & 4){
                                if(f->obj->smooth) colrgb = colorToRGB(lerp(col, col2, col3, double(w1[i])/w0[i], double(w2[i])/w0[i], double(w3[i])/w0[i]));
                                else if(colrgb>>24) colrgb = colorToRGB(calcLighting(fcenter, f->normal, *f->obj->mat, sc));
                                imgrasters->colbuffer[pint.y*w + pint.x+i] = colrgb;
                            }
                            imgrasters->zbuffer[w*pint.y+pint.x+i] = z[i];
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

void zBufferDraw(raster *imgrasters, scene *sc){
#ifdef USE_VECTOR
    generate_maps_vector(5, imgrasters, sc);
#else
    generate_maps(5, imgrasters, sc);
#endif
}

void paintNormalMap(raster *imgrasters, scene *sc){
    generate_maps_vector(3, imgrasters, sc);
    std::copy(imgrasters->normbuffer, imgrasters->normbuffer+(imgrasters->width()*imgrasters->height()), imgrasters->colbuffer);
}

//TODO not really SSAO, should probably fix some tweaks
void SSAO(raster *imgrasters, scene *sc)
{
    int w=imgrasters->width(), h=imgrasters->height();
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
            // * color(1,1,1) *
            imgrasters->colbuffer[y*w + x] = colorToRGB(RGBToColor(imgrasters->colbuffer[y*w+x])*clamp(ao*2, 0, 1));
        }
    }
    
}

// raytraces an image, breaking it up into tiles of dimension tilesize x tilesize
void rayTraceUnthreaded(raster *imgrasters, scene *sc, int tilesize){
    num_rays_traced = 0;
    int i, j, x, y, xmax, ymax, w=imgrasters->width(), h=imgrasters->height();
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
                    col = traceRay(r, sc);
                    colrgb = colorToRGB(col);
                    imgrasters->colbuffer[y*w + x] = colrgb;
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

// path traces an image
void pathTraceUnthreaded(raster *imgrasters, scene *sc, int tilesize){
    num_rays_traced = 0;
    clock_t begin = clock();
    
    int i, j, x, y, xmax, ymax, w=imgrasters->width(), h=imgrasters->height(), s;
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
                    for(s=1, totalcol=color(); s<=PATH_TRACE_SAMPLES; s++){
                        r = sc->cam->castRay(x, y, w, h);
                        totalcol += tracePath(r, sc);
                    }
                    colrgb = colorToRGB(totalcol*(1.0/PATH_TRACE_SAMPLES));
                    imgrasters->colbuffer[y*w + x] = colrgb;
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

color rayTracePixel(int x, int y, int w, int h, scene *sc){
    return traceRay(sc->cam->castRay(x, y, w, h), sc);
}

color pathTracePixel(int x, int y, int w, int h, scene *sc){
    int s;
    color totalcol;
    for(s=1, totalcol=color(); s<=PATH_TRACE_SAMPLES; s++){
        totalcol += tracePath(sc->cam->castRay(x, y, w, h), sc);
    }
    return totalcol*(1.0/s);
}

color ambOccPixel(int x, int y, int w, int h, scene *sc){
    ray r = sc->cam->castRay(x, y, w, h);
    double ao = ambientOcclusion(r, sc);
    ao = clamp(ao*2, 0, 1);
    return color(ao, ao, ao);
}

color ray_march_pixel(int x, int y, int w, int h, scene *sc){
    return rayTraceDistanceField(sc->cam->castRay(x, y, w, h), sc, 50);
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