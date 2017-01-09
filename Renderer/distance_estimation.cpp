//
//  fractal.cpp
//  Renderer
//
//  Created by Raphael Kargon on 1/5/17.
//  Copyright © 2017 Raphael Kargon. All rights reserved.
//

#include "distance_estimation.h"

// Smoothing is done with a polynomial function that has a continuous derivative.
double smooth_min(double a, double b, double k){
    double h = clamp(0.5 + 0.5 * (b-a)/k, 0, 1);
    return lerp(b, a, h) - k*h*(1.0-h);
}

double smooth_max(double a, double b, double k){
    return -smooth_min(-a, -b, k);
}

/* Distance Estimator Primitives */

distance_estimator de_plane(vertex normal){
    return [=](const vertex& v) -> double {return dot(normal, v);};
}

distance_estimator de_sphere(vertex center, double radius){
    return [=](const vertex& v) -> double {return (v-center).len() - radius;};
}

//distance_estimator de_box(const bounds& bounding_box){
//    return [=](const vertex& v) -> double {return (v-center).len() - radius;};
//}

distance_estimator de_torus(double r1, double r2){
    return [=](const vertex& v) -> double {
        // project vertex onto ring
        vertex ring_pt = (vertex{1,1,0} * v).unitvect() * r1;
        return (v - ring_pt).len() - r2;
    };
}

distance_estimator de_sierpinski_tetrahedron(int num_iterations){
    return [=](const vertex& v) -> double {
        vertex z = v;
        vertex offset = {1,1,1};
        int n;
        for(n=0; n<num_iterations; n++){
            // Fold point across planes to make sure it's in the correct (top) octant of tetrahedron
            if (z.x + z.y < 0) {
                z.set(-z.y, -z.x, z.z);
            }
            if (z.y + z.z < 0) {
                z.set(z.x, -z.z, -z.y);
            }
            if (z.z + z.x < 0) {
                z.set(-z.z, z.y, -z.x);
            }
            // rotate point
            //z = 2 * rotate(z-offset, {0,0,1}, 0.1) + offset;
            
            // Scale points up
            z = 2.0*z - offset;
        }
        double p = pow(2.0, -num_iterations);
        return z.len() * p;
    };
}

distance_estimator de_menger(double scale, vertex center, int num_iterations){
    return [=](const vertex& v) -> double {
        vertex z = v;
        int n;
        for (n=0; n<num_iterations; n++){
            // rotate 1
//            z = rotate(z, {0,-1,1}, 0.2);
            
            z = abs(z);
            if (z.x - z.y < 0){
                z.set(z.y, z.x, z.z);
            }
            if (z.x - z.z < 0){
                z.set(z.z, z.y, z.x);
            }
            if (z.y - z.z < 0){
                z.set(z.x, z.z, z.y);
            }
            
            // rotate 2
//            z = rotate(z, {0,0,1}, 0.1);
            
            // TODO fix scaling
            double tmp_z = z.z;
            z = scale*(z - center) + center;
            z.z = tmp_z * scale;
            if(z.z > 0.5*center.z*(scale-1)){
                z.z -= center.z*(scale-1);
            }
            
        }
        return z.len() * pow(scale, -n);
    };
}

//Menger3(x,y,z){
//    r=x*x+y*y+z*z;
//    for(i=0;i<MI && r<bailout;i++){
//        rotate1(x,y,z);
//        
//        x=abs(x);y=abs(y);z=abs(z);
//        if(x-y<0){x1=y;y=x;x=x1;}
//        if(x-z<0){x1=z;z=x;x=x1;}
//        if(y-z<0){y1=z;z=y;y=y1;}
//        
//        rotate2(x,y,z);
//        
//        x=scale*x-CX*(scale-1);
//        y=scale*y-CY*(scale-1);
//        z=scale*z;
//        if(z>0.5*CZ*(scale-1)) z-=CZ*(scale-1);
//        
//        r=x*x+y*y+z*z;
//    }
//    return (sqrt(x*x+y*y+z*z)-2)*scale^(-i);
//}

/* Distance Estimator Operations */

distance_estimator de_union(distance_estimator de1, distance_estimator de2){
    return [=](const vertex& v) -> double {return std::min(de1(v), de2(v));};
}

distance_estimator operator||(distance_estimator de1, distance_estimator de2){
    return de_union(de1, de2);
}

distance_estimator de_blend(distance_estimator de1, distance_estimator de2, double k){
    return [=](const vertex& v) -> double {return smooth_min(de1(v), de2(v), k);};
}

distance_estimator de_intersect(distance_estimator de1, distance_estimator de2){
    return [=](const vertex& v) -> double {return std::max(de1(v), de2(v));};
}

distance_estimator operator&&(distance_estimator de1, distance_estimator de2){
    return de_intersect(de1, de2);
}

distance_estimator de_intersect_blend(distance_estimator de1, distance_estimator de2, double k){
    return [=](const vertex& v) -> double {return smooth_max(de1(v), de2(v), k);};
}

distance_estimator de_subtract(distance_estimator de1, distance_estimator de2){
    return [=](const vertex& v) -> double {return std::max(de1(v), -de2(v));};
}

distance_estimator de_twist(distance_estimator de){
    return [=](const vertex& v) -> double {
        double c = cos(4*v.y);
        double s = sin(4*v.y);
        return de({c*v.x - s*v.z, s*v.x + c*v.z, v.y});
    };
}

vertex estimate_normal(const vertex& v, const distance_estimator& obj, double epsilon){
    vertex dx{epsilon, 0, 0};
    vertex dy{0, epsilon, 0};
    vertex dz{0, 0, epsilon};
    vertex normal = vertex(obj(v+dx) - obj(v-dx),
                           obj(v+dy) - obj(v-dy),
                           obj(v+dz) - obj(v-dz)).unitvect();
    return normal;
}

bool ray_march(const ray& r, const distance_estimator& obj, double *final_dist, int *steps, int max_ray_steps, double epsilon){
    double t = epsilon*1.5;
    for(int s=0; s<max_ray_steps; s++){
        vertex v = r.org + t*r.dir;
        double distance = obj(v);
        t += distance;
        if (distance <= epsilon){
            if (steps != nullptr){
                *steps = s;
            }
            if (final_dist != nullptr){
                *final_dist = t;
            }
            return (t > epsilon);
        }
    }
    return false;
}
