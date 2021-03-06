//
//  rendering.cpp - the magic happens here!
//  Renderer
//
//  Created by Raphael Kargon on 6/8/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "rendering.h"

#include "sampling.h"

#include <glm/gtx/io.hpp>
#include <glm/gtx/norm.hpp>

#include <random>

int num_rays_traced = 0;

// PRECONDITION: N is normalized
// Calculate lighting for triangle rasterization.
// Given a vertex, a normal, and a material, will calculate the color of it,
// taking into account incident light from lamps.
// Does not take into account other geometry which could cause shadows,
// reflections etc.
color calc_lighting(const vertex &v, const vertex &n, const material &mat,
                    const scene &sc) {
  color lightcol(0), speclightcol(0);
  for (const lamp &l : sc.lamps) {
    vertex lampvect = l.loc - v;
    vertex lampvnorm = glm::normalize(lampvect);
    double dotprod = glm::dot(n, lampvnorm);
    vertex view = glm::normalize(sc.cam.view_vector(v));
    // if lamp and view are on different sides of face, then one is looking at
    // underside of face)
    if (dotprod * glm::dot(n, view) > 0) {
      continue;
    }
    double dstsqr = glm::length2(lampvect);
    if (dstsqr == 0) {
      for (int c = 0; c < 3; ++c) {
        if (l.col[c]) {
          lightcol[c]++;
        }
      }
      continue;
    }
    // Phong specular reflection
    vertex refl = glm::reflect(lampvnorm, n);
    double spec_intensity =
        l.intensity *
        std::fmax(0, std::pow(glm::dot(view, refl), mat.spec_hardness));
    double diff_intensity = l.intensity * std::fabs(dotprod);
    // calc falloff
    diff_intensity /= dstsqr;
    spec_intensity /= dstsqr;
    lightcol += l.col * diff_intensity;
    speclightcol += l.col * spec_intensity;
  }
  color col = mat.diff_col * lightcol;
  color spcol = mat.spec_col * speclightcol;
  return col + spcol;
}

// Traces a ray from the camera through a scene, using ray-tracing, up to a
// certain depth
color trace_ray(const ray &viewray, const scene &sc, unsigned int depth,
                unsigned int max_depth) {
  num_rays_traced++;
  vertex tuv;
  const face *f = kdtree::ray_tree_intersect(&sc.kdt, viewray, false, &tuv);
  if (f == nullptr) {
    return sc.w->get_color(viewray);
  } else {
    const mesh &obj = *f->obj;
    const material &mat = sc.materials[obj.mat_id];
    vertex v = viewray.org + viewray.dir * tuv[0]; // calculate vertex location

    // interpolate texture coordinates
    double txu, txv;
    // XXX TODO
    // txu = f->vertices[0]->tex_u * (1 - tuv[1] - tuv[2]) +
    //       f->vertices[1]->tex_u * tuv[1] + f->vertices[2]->tex_u * tuv[2];
    // txv = f->vertices[0]->tex_v * (1 - tuv[1] - tuv[2]) +
    //       f->vertices[1]->tex_v * tuv[1] + f->vertices[2]->tex_v * tuv[2];
    txu = txv = 0;

    // calculate normal
    vertex n;
    if (obj.smooth) {
      n = (obj.vertex_normals[f->v[0]]) * (1 - tuv[1] - tuv[2]) +
          (obj.vertex_normals[f->v[1]]) * tuv[1] +
          (obj.vertex_normals[f->v[2]]) * tuv[2];
      normalize_in_place(n);
    } else {
      n = f->normal;
    }
    // n = sc.obj.mat->get_normal(txu, txv);

    double ndotray = dot(n, viewray.dir);

    /* SPECULAR & DIFFUSE LIGHTING */
    color lightcol(0), speclightcol(0);
    for (const lamp &l : sc.lamps) {
      vertex lampvect = l.loc - v;
      vertex lampvnorm = glm::normalize(lampvect);
      // make sure lamp is on same side of face as view
      if (glm::dot(n, lampvect) * ndotray > 0) {
        continue;
      }
      ray lampray(v, lampvnorm);
      if (kdtree::ray_tree_intersect(&sc.kdt, lampray, true, &tuv) != nullptr) {
        continue;
      }
      double dotprod = glm::dot(lampvnorm, n);
      double dstsqr = glm::length2(lampvect);
      if (dstsqr == 0) {
        // if lamp is on the face's center, add full brightness to each color
        // that the lamp emits
        for (int c = 0; c < 3; ++c) {
          if (l.col[c]) {
            lightcol[c]++;
          }
        }
        continue;
      }
      // Phong shading
      vertex lamprefl = glm::reflect(lampvnorm, n);
      double spec_intensity =
          l.intensity * std::fmax(0, std::pow(glm::dot(viewray.dir, lamprefl),
                                              mat.spec_hardness));
      double diff_intensity = l.intensity * std::fabs(dotprod);
      // calculate falloff
      switch (l.falloff) {
      case 2:
        diff_intensity /= dstsqr;
        spec_intensity /= dstsqr;
        break;
      case 1:
        diff_intensity /= std::sqrt(dstsqr);
        spec_intensity /= std::sqrt(dstsqr);
        break;
      case 0:
        break;
      }
      lightcol += l.col * diff_intensity;
      speclightcol += l.col * spec_intensity;
    }

    color totcol = mat.get_color(txu, txv) * lightcol +
                   mat.get_spec_col(txu, txv) * speclightcol;

    /* Ambient lighting from sky */
    if (sc.w->ambient_intensity > 0) {
      ray normal_ray{v, n};
      totcol += sc.w->ambient_intensity * sc.w->get_color(normal_ray);
    }

    /* REFLECTION & REFRACTION */
    // an approximation. Assumes normals point outside and doesn't really deal
    // with with concentric/intersecting objects.Currently assumes 'outside' of
    // every object is air.
    // also doesn't do fresnel formula, reflection and refraction are handled
    // separately, except for total internal reflection
    if (depth < max_depth) {
      vertex refl = glm::reflect(viewray.dir, n);
      if (mat.alpha < 1) {
        double n1, n2;
        vertex transray;
        if (ndotray < 0) {
          n1 = 1;
          n2 = mat.ior;
        } else {
          n1 = mat.ior;
          n2 = 1;
        }
        vertex raytang = viewray.dir - n * ndotray;
        vertex transtang = raytang * (n1 / n2);
        double transsinsquared = glm::length2(transtang);
        if (transsinsquared > 1)
          transray = refl; // total reflection
        else {
          vertex transnorm = n * signum(ndotray) * sqrt(1 - transsinsquared);
          transray = transnorm + transtang;
        }
        color transcol = trace_ray(ray(v, transray), sc, depth + 1, max_depth);
        totcol = glm::mix(transcol, totcol, mat.alpha);
      }
      if (mat.refl_intensity > 0) {
        color refcol =
            mat.spec_col * trace_ray(ray(v, refl), sc, depth + 1, max_depth);
        totcol = glm::mix(totcol, refcol, mat.refl_intensity);
      }
    }
    return totcol;
  }
}

// TODO transparent objects cast shadows

// Traces a path through a scene, up to the given depth.
// Used for path-tracing, and calculates indirect lighting.
color trace_path(const ray &viewray, const scene &sc, unsigned int depth,
                 unsigned int max_depth) {
  if (depth > max_depth) {
    return color(0);
  }

  num_rays_traced++;
  vertex tuv;
  const face *f = kdtree::ray_tree_intersect(&sc.kdt, viewray, false, &tuv);
  if (f == nullptr) {
    return sc.w->get_color(viewray);
  }
  const mesh &obj = *f->obj;
  vertex v = viewray.org + viewray.dir * tuv[0]; // calculate vertex location

  // calculate normal
  vertex n;
  if (obj.smooth) {
    n = (obj.vertex_normals[f->v[0]]) * (1 - tuv[1] - tuv[2]) +
        (obj.vertex_normals[f->v[1]]) * tuv[1] +
        (obj.vertex_normals[f->v[2]]) * tuv[2];
    normalize_in_place(n);
  } else {
    n = f->normal;
  }

  auto emitter = dynamic_cast<const EmissionBSDF *>(obj.bsdf);
  if (emitter) {
    return emitter->col;
  }
  double pdf_rr = 0.9; // TODO smarter russian roullette
  if (std::uniform_real_distribution<float>(0, 1)(generator) > pdf_rr) {
    return color(0);
  }

  auto wi_pdf = obj.bsdf->sample_direction(n, -viewray.dir);
  vertex wi = wi_pdf.first;
  double pdf = wi_pdf.second;
  color incoming_light = trace_path(ray(v + EPSILON * wi, wi), sc, depth + 1);
  return incoming_light * obj.bsdf->bsdf(-wi, -viewray.dir) * glm::dot(wi, n) /
         (pdf * pdf_rr);

  // // calculate incident light
  // vertex inc_dir = obj.bsdf->getIncidentDirection(n, viewray.dir);
  // color inc_col(0);
  // // TODO use russian roullette
  // if (!dynamic_cast<const EmissionBSDF *>(obj.bsdf)) {
  //   inc_col = trace_path(ray(v, inc_dir), sc, depth + 1);
  // }
  // // calculate returned light
  // color return_col = obj.bsdf->getLight(inc_col, inc_dir, n, viewray.dir);
  // return return_col;
}

// Calculates ambient occlusion for a ray.
double ambient_occlusion(const ray &viewray, const scene &sc, int samples) {
  vertex tuv;
  const face *f = kdtree::ray_tree_intersect(&sc.kdt, viewray, false, &tuv);
  if (f == nullptr)
    return 1;
  else {
    const mesh &obj = *f->obj;
    vertex v = viewray.org + viewray.dir * tuv[0]; // calculate vertex location
    // calculate normal
    vertex n;
    if (obj.smooth) {
      n = (obj.vertex_normals[f->v[0]]) * (1 - tuv[1] - tuv[2]) +
          (obj.vertex_normals[f->v[1]]) * tuv[1] +
          (obj.vertex_normals[f->v[2]]) * tuv[2];
      normalize_in_place(n);
    } else {
      n = f->normal;
    }
    // double ndotray = dot(n, viewray.dir);

    double occ_amount = 0; // amount of ambient occlusion, ie how much current
                           // point is illuminated by the background.
    ray testray(v, vertex(0, 0, 0));
    for (int i = 1; i <= samples; i++) {
      testray.dir = random_direction();
      if (dot(testray.dir, n) < 0)
        continue;
      else if (kdtree::ray_tree_intersect(&sc.kdt, testray, true, nullptr) ==
               nullptr)
        occ_amount++;
    }
    occ_amount /= samples;
    return occ_amount;
  }
}

color raytrace_distance_field(const ray &viewray, const scene &sc,
                              int num_steps, unsigned int depth,
                              unsigned int max_depth) {
  num_rays_traced++;
  double t;
  int steps;
  if (!ray_march(viewray, sc.de_obj, &t, &steps, num_steps)) {
    return sc.w->get_color(viewray);
  } else {
    vertex v = viewray.org + t * viewray.dir;
    double ambient_occlusion = 1.0 - (steps / (double)num_steps);
    vertex n = estimate_normal(v, sc.de_obj);
    double ndotray = dot(n, viewray.dir);

    /* SPECULAR & DIFFUSE LIGHTING */
    color lightcol, speclightcol;
    for (const lamp &l : sc.lamps) {
      vertex lampvect = l.loc - v;
      vertex lampvnorm = glm::normalize(lampvect);
      if (dot(n, lampvect) * ndotray > 0)
        continue; // make sure lamp is on same side of face as view
      ray lampray(v, lampvnorm);
      if (ray_march(lampray, sc.de_obj, nullptr, nullptr, num_steps)) {
        continue;
      }
      double dotprod = glm::dot(lampvnorm, n);
      double dstsqr = glm::length2(lampvect);
      if (dstsqr == 0) {
        // if lamp is on the face's center, add full brightness to each color
        // that the lamp emits
        if (l.col.r)
          lightcol.r++;
        if (l.col.g)
          lightcol.g++;
        if (l.col.b)
          lightcol.b++;
        continue;
      }
      // Phong shading
      vertex lamprefl = glm::reflect(lampvnorm, n);
      double spec_intensity =
          l.intensity *
          fmax(0, pow(dot(viewray.dir, lamprefl), sc.de_mat.spec_hardness));
      double diff_intensity = l.intensity * fabs(dotprod);
      // calculate falloff
      switch (l.falloff) {
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
      lightcol += l.col * diff_intensity;
      speclightcol += l.col * spec_intensity;
    }
    color totcol =
        sc.de_mat.diff_col * lightcol + sc.de_mat.spec_col * speclightcol;

    /* Ambient lighting from sky */
    if (sc.w->ambient_intensity > 0) {
      ray normal_ray{v, n};
      totcol += sc.w->ambient_intensity * sc.w->get_color(normal_ray) *
                ambient_occlusion;
    }

    /* REFLECTION & REFRACTION */
    // An approximation. Assumes normals point outside and doesn't really deal
    // with with concentric/intersecting objects.
    // Currently assumes 'outside' of every object is air.
    // also doesn't do fresnel formula, reflection and refraction are handled
    // separately, except for total internal reflection
    if (depth < max_depth) {
      vertex refl = glm::reflect(viewray.dir, n);
      if (sc.de_mat.alpha < 1) {
        double n1, n2;
        vertex transray;
        if (ndotray < 0) {
          n1 = 1;
          n2 = sc.de_mat.ior;
        } else {
          n1 = sc.de_mat.ior;
          n2 = 1;
        }
        vertex raynorm = n * ndotray;
        vertex raytang = viewray.dir - raynorm;
        vertex transtang = raytang * (n1 / n2);
        double transsinsquared = glm::length2(transtang);
        if (transsinsquared > 1) {
          transray = refl;
        } // total reflection
        else {
          vertex transnorm = n * signum(ndotray) * sqrt(1 - transsinsquared);
          transray = transnorm + transtang;
        }
        // TODO
        color transcol = trace_ray(ray(v, transray), sc, depth + 1);
        totcol = glm::mix(transcol, totcol, sc.de_mat.alpha);
      }
      if (sc.de_mat.refl_intensity > 0) {
        color refcol =
            sc.de_mat.spec_col * trace_ray(ray(v, refl), sc, depth + 1);
        totcol = glm::mix(totcol, refcol, sc.de_mat.refl_intensity);
      }
    }
    return totcol;
  }
}

/* Rasterization */

// generates color, normal, and depth maps
// mapflags bit flags:
// 1 - depth map - (will always be generated anyway, needed for other maps)
// 2 - normal map
// 4 - color map
// reference implementation, no vector operations
void generate_maps(int mapflags, raster &imgrasters, const scene &sc) {
  int w = imgrasters.width(), h = imgrasters.height();
  double z, z1, z2, z3, dz21, dz31;
  int minx, miny, maxx, maxy;
  int A12, A23, A31, B12, B23, B31;
  int w0, w1, w2, w3, w1_row, w2_row, w3_row;
  int wsgn;
  vertex fcenter;
  point_2d<double> p1, p2, p3, p;
  point_2d<int> p1int, p2int, p3int, pint;
  color col, col2, col3;
  vertex norm, norm2, norm3, normtmp;
  uint colrgb = 1;

  std::fill_n(imgrasters.zbuffer.get(), imgrasters.size(), 1);
  if (mapflags & 2)
    std::fill_n(imgrasters.normbuffer.get(), imgrasters.size(), 0xffffff);
  if (mapflags & 4)
    std::fill_n(imgrasters.colbuffer.get(), imgrasters.size(), 0xffffff);

  for (const auto &obj_ptr : sc.objects) {
    for (const face &f : obj_ptr->faces) {
      // TODO per-face materials
      const material &mat = sc.materials[obj_ptr->mat_id];

      // get pixels of vertices
      p1 = sc.cam.project_vertex(f.get_vert(0), w, h);
      p2 = sc.cam.project_vertex(f.get_vert(1), w, h);
      p3 = sc.cam.project_vertex(f.get_vert(2), w, h);
      if (isnan(p1.x) || isnan(p1.y) || isnan(p2.x) || isnan(p2.y) ||
          isnan(p3.x) || isnan(p3.y))
        continue;
      p1int.x = (int)p1.x;
      p1int.y = (int)p1.y;
      p2int.x = (int)p2.x;
      p2int.y = (int)p2.y;
      p3int.x = (int)p3.x;
      p3int.y = (int)p3.y;

      // z values = (z-min)/(max-min)
      z1 = (sc.cam.vertex_depth(f.get_vert(0)) - sc.cam.mindist) /
           (sc.cam.maxdist - sc.cam.mindist);
      z2 = (sc.cam.vertex_depth(f.get_vert(1)) - sc.cam.mindist) /
           (sc.cam.maxdist - sc.cam.mindist);
      z3 = (sc.cam.vertex_depth(f.get_vert(2)) - sc.cam.mindist) /
           (sc.cam.maxdist - sc.cam.mindist);

      // store difference values. Makes interpolation later on slightly faster
      dz21 = z2 - z1;
      dz31 = z3 - z1;

      if (mapflags & 4) {
        fcenter = f.center();
        colrgb = 1 << 24; // unitialized color, largest byte is non-zero
      }

      // smooth shading, get vertex colors
      if ((mapflags & (4 + 2)) && f.obj->smooth) {
        norm = obj_ptr->vertex_normals[f.v[0]];
        norm2 = obj_ptr->vertex_normals[f.v[1]];
        norm3 = obj_ptr->vertex_normals[f.v[2]];
        col = calc_lighting(f.get_vert(0), norm, mat, sc);
        col2 = calc_lighting(f.get_vert(1), norm2, mat, sc);
        col3 = calc_lighting(f.get_vert(2), norm3, mat, sc);
      }

      // triangle bounding box
      minx = std::min(std::min(p1int.x, p2int.x), p3int.x);
      miny = std::min(std::min(p1int.y, p2int.y), p3int.y);
      maxx = std::max(std::max(p1int.x, p2int.x), p3int.x);
      maxy = std::max(std::max(p1int.y, p2int.y), p3int.y);
      if (minx > w - 1 || maxx < 0 || miny > h - 1 || maxy < 0)
        continue; // face is off screen

      // clipint to screen
      minx = std::max(minx, 0);
      maxx = std::min(maxx, w - 1);
      miny = std::max(miny, 0);
      maxy = std::min(maxy, h - 1);

      // triangle edge setup
      A12 = p1int.y - p2int.y;
      A23 = p2int.y - p3int.y;
      A31 = p3int.y - p1int.y;
      B12 = p2int.x - p1int.x;
      B23 = p3int.x - p2int.x;
      B31 = p1int.x - p3int.x;

      // initial barycentric coordinates at corner
      pint = point_2d<int>(minx, miny);
      w1_row = orient_2d(p2int, p3int, pint);
      w2_row = orient_2d(p3int, p1int, pint);
      w3_row = orient_2d(p1int, p2int, pint);
      w0 = orient_2d(p1int, p2int, p3int);
      if (w0 == 0) {
        continue;
      }
      wsgn = signum(w0);

      // rasterize
      for (pint.y = miny; pint.y <= maxy; pint.y++) {
        w1 = w1_row;
        w2 = w2_row;
        w3 = w3_row;
        for (pint.x = minx; pint.x <= maxx; pint.x++) {
          if ((signum(w1) == wsgn || !w1) && (signum(w2) == wsgn || !w2) &&
              (signum(w3) == wsgn || !w3)) {
            // interpolate z value
            z = z1 + w2 * dz21 / w0 + w3 * dz31 / w0;
            if (z < imgrasters.zbuffer[w * pint.y + pint.x]) {
              if (mapflags & 2) {
                if (obj_ptr->smooth) {
                  normtmp = lerp(norm, norm2, norm3, double(w1) / w0,
                                 double(w2) / w0, double(w3) / w0);
                } else {
                  normtmp = f.normal;
                }
                imgrasters.normbuffer[pint.y * w + pint.x] =
                    normal_to_rgb(normtmp);
              }
              if (mapflags & 4) {
                if (f.obj->smooth) {
                  colrgb = color_to_rgb(lerp(col, col2, col3, double(w1) / w0,
                                             double(w2) / w0, double(w3) / w0));
                } else if (colrgb >> 24) {
                  // TODO per-face materials
                  colrgb =
                      color_to_rgb(calc_lighting(fcenter, f.normal, mat, sc));
                }
                imgrasters.colbuffer[pint.y * w + pint.x] = colrgb;
              }
              imgrasters.zbuffer[w * pint.y + pint.x] = z;
            }
          }
          w1 += A23;
          w2 += A31;
          w3 += A12;
        }
        w1_row += B23;
        w2_row += B31;
        w3_row += B12;
      }
    }
  }
}

// generates color, normal, and depth maps
// mapflags bit flags:
// 1 - depth map - (will always be generated anyway, needed for other maps)
// 2 - normal map
// 4 - color map
void generate_maps_vector(int mapflags, raster &imgrasters, const scene &sc) {
  int i;
  int w = imgrasters.width(), h = imgrasters.height();
  int minx, miny, maxx, maxy;
  __v4sf z, z1, z2, z3, dz21, dz31; // interpolated z values
  __v4si w0, w1, w2, w3, w1_row, w2_row, w3_row, wsgn;
  __v4si pxmask; // whether each pixel is inside the triangle
  vertex fcenter;
  point_2d<double> p1, p2, p3;
  point_2d<int> p1int, p2int, p3int, pint;
  color col, col2, col3;
  vertex norm, norm2, norm3, normtmp;
  uint colrgb = 1;

  std::fill_n(imgrasters.zbuffer.get(), imgrasters.size(), 1);
  if (mapflags & 2) {
    std::fill_n(imgrasters.normbuffer.get(), imgrasters.size(), 0xffffff);
  }
  if (mapflags & 4) {
    std::fill_n(imgrasters.colbuffer.get(), imgrasters.size(), 0xffffff);
  }

  for (const auto &obj_ptr : sc.objects) {
    for (const face &f : obj_ptr->faces) {
      // TODO per-face materials
      const material &mat = sc.materials[obj_ptr->mat_id];
      // get pixels of vertices
      p1 = sc.cam.project_vertex(f.get_vert(0), w, h);
      p2 = sc.cam.project_vertex(f.get_vert(1), w, h);
      p3 = sc.cam.project_vertex(f.get_vert(2), w, h);

      if (isnan(p1.x) || isnan(p1.y) || isnan(p2.x) || isnan(p2.y) ||
          isnan(p3.x) || isnan(p3.y)) {
        continue;
      }

      p1int.x = (int)p1.x;
      p1int.y = (int)p1.y;
      p2int.x = (int)p2.x;
      p2int.y = (int)p2.y;
      p3int.x = (int)p3.x;
      p3int.y = (int)p3.y;

      // z values
      z1 = _mm_set1_ps((sc.cam.vertex_depth(f.get_vert(0)) - sc.cam.mindist) /
                       (sc.cam.maxdist - sc.cam.mindist));
      z2 = _mm_set1_ps((sc.cam.vertex_depth(f.get_vert(1)) - sc.cam.mindist) /
                       (sc.cam.maxdist - sc.cam.mindist));
      z3 = _mm_set1_ps((sc.cam.vertex_depth(f.get_vert(2)) - sc.cam.mindist) /
                       (sc.cam.maxdist - sc.cam.mindist));

      // store difference values. Makes interpolation later on slightly faster
      dz21 = z2 - z1;
      dz31 = z3 - z1;

      if (mapflags & 4) {
        fcenter = f.center();
        colrgb = 1 << 24; // unitialized color, largest byte is non-zero
      }

      // smooth shading, get vertex colors
      if ((mapflags & (4 + 2)) && obj_ptr->smooth) {
        norm = obj_ptr->vertex_normals[f.v[0]];
        norm2 = obj_ptr->vertex_normals[f.v[1]];
        norm3 = obj_ptr->vertex_normals[f.v[2]];
        col = calc_lighting(f.get_vert(0), norm, mat, sc);
        col2 = calc_lighting(f.get_vert(1), norm2, mat, sc);
        col3 = calc_lighting(f.get_vert(2), norm3, mat, sc);
      }

      // triangle bounding box
      minx = std::min(std::min(p1int.x, p2int.x), p3int.x);
      miny = std::min(std::min(p1int.y, p2int.y), p3int.y);
      maxx = std::max(std::max(p1int.x, p2int.x), p3int.x);
      maxy = std::max(std::max(p1int.y, p2int.y), p3int.y);
      if (minx > w - 1 || maxx < 0 || miny > h - 1 || maxy < 0)
        continue; // face is off screen

      // clip to screen
      minx = std::max(minx, 0);
      maxx = std::min(maxx, w - 1);
      miny = std::max(miny, 0);
      maxy = std::min(maxy, h - 1);

      // triangle edge setup
      pint = point_2d<int>(minx, miny);
      edge_vect e12, e23, e31;

      // initial barycentric coordinates at corner
      w1_row = e23.init(p2int, p3int, pint);
      w2_row = e31.init(p3int, p1int, pint);
      w3_row = e12.init(p1int, p2int, pint);
      w0 = _mm_set1_epi32(orient_2d(p1int, p2int, p3int));
      // if w0 is zero continue
      if (_mm_movemask_epi8(_mm_cmpeq_epi32(w0, zeroveci)) > 0) {
        continue;
      }
      wsgn = _mm_cmpgt_epi32(w0, zeroveci);

      // rasterize
      for (pint.y = miny; pint.y <= maxy; pint.y += edge_vect::step_y) {
        w1 = w1_row;
        w2 = w2_row;
        w3 = w3_row;
        for (pint.x = minx; pint.x <= maxx; pint.x += edge_vect::step_x) {
          // each item in pxmask vector is >0 iff corresponding pixel should be
          // drawn
          // check if w1,2,3 have the same sign as w0, or are 0 themselves
          pxmask =
              _mm_or_si128(_mm_cmpeq_epi32(_mm_cmpgt_epi32(w1, zeroveci), wsgn),
                           _mm_cmpeq_epi32(w1, zeroveci));
          pxmask = _mm_and_si128(
              pxmask,
              _mm_or_si128(_mm_cmpeq_epi32(_mm_cmpgt_epi32(w2, zeroveci), wsgn),
                           _mm_cmpeq_epi32(w2, zeroveci)));
          pxmask = _mm_and_si128(
              pxmask,
              _mm_or_si128(_mm_cmpeq_epi32(_mm_cmpgt_epi32(w3, zeroveci), wsgn),
                           _mm_cmpeq_epi32(w3, zeroveci)));
          // interpolate z value
          // z = z1 + w2*dz21/w0 + w3*dz31/w0
          z = z1 + _mm_div_ps(_mm_mul_ps(_mm_cvtepi32_ps(w2), dz21),
                              _mm_cvtepi32_ps(w0)) +
              _mm_div_ps(_mm_mul_ps(_mm_cvtepi32_ps(w3), dz31),
                         _mm_cvtepi32_ps(w0));

          for (i = 0; i < 4; i++) {
            if (pint.x + i < w && pxmask[i] != 0 &&
                z[i] < imgrasters.zbuffer[w * pint.y + pint.x + i]) {
              if (mapflags & 2) {
                if (obj_ptr->smooth) {
                  normtmp = lerp(norm, norm2, norm3, double(w1[i]) / w0[i],
                                 double(w2[i]) / w0[i], double(w3[i]) / w0[i]);
                } else {
                  normtmp = f.normal;
                }
                imgrasters.normbuffer[pint.y * w + pint.x + i] =
                    normal_to_rgb(normtmp);
              }
              if (mapflags & 4) {
                if (obj_ptr->smooth) {
                  colrgb = color_to_rgb(
                      lerp(col, col2, col3, double(w1[i]) / w0[i],
                           double(w2[i]) / w0[i], double(w3[i]) / w0[i]));
                } else if (colrgb >> 24) {
                  colrgb =
                      color_to_rgb(calc_lighting(fcenter, f.normal, mat, sc));
                }
                imgrasters.colbuffer[pint.y * w + pint.x + i] = colrgb;
              }
              imgrasters.zbuffer[w * pint.y + pint.x + i] = z[i];
            }
          }

          w1 += e23.one_step_x;
          w2 += e31.one_step_x;
          w3 += e12.one_step_x;
        }
        w1_row += e23.one_step_y;
        w2_row += e31.one_step_y;
        w3_row += e12.one_step_y;
      }
    }
  }
}

void zbuffer_draw(raster &imgrasters, const scene &sc) {
#ifdef USE_VECTOR
  generate_maps_vector(5, imgrasters, sc);
#else
  generate_maps(5, imgrasters, sc);
#endif
}

void paint_normal_map(raster &imgrasters, const scene &sc) {
  generate_maps_vector(3, imgrasters, sc);
  std::copy(imgrasters.normbuffer.get(),
            imgrasters.normbuffer.get() +
                (imgrasters.width() * imgrasters.height()),
            imgrasters.colbuffer.get());
}

// TODO not really SSAO, should probably fix some tweaks
void SSAO(raster &imgrasters, const scene &sc) {
  int w = imgrasters.width(), h = imgrasters.height();
  generate_maps_vector(7, imgrasters, sc); // generate depth and normal maps
  double ao, z, ztmp;
  vertex v, dv, vtmp, n;
  ray r, rtmp;
  int x, y, dy, dx, dir, nsamples;

  for (y = 0; y < h; y++) {
    for (x = 0; x < w; x++) {
      z = imgrasters.zbuffer[y * w + x];
      if (z == 1) {
        continue;
      }
      z = z * (sc.cam.maxdist - sc.cam.mindist) + sc.cam.mindist;
      r = sc.cam.cast_ray(x, y, w, h);
      v = r.org + r.dir * z;
      n = rgb_to_normal(imgrasters.normbuffer[y * w + x]);

      ao = 0;
      nsamples = 0;
      for (dir = 0; dir < 20; dir++) {
        dx = rand() % 10 - 5;
        dy = rand() % 10 - 5;
        if (x + dx < 0 || x + dx >= w) {
          continue;
        }
        if (y + dy < 0 || y + dy >= h) {
          continue;
        }

        ztmp = imgrasters.zbuffer[(y + dy) * w + (x + dx)];
        if (ztmp == 1) {
          continue;
        }
        ztmp = ztmp * (sc.cam.maxdist - sc.cam.mindist) + sc.cam.mindist;
        rtmp = sc.cam.cast_ray(x, y, w, h);
        vtmp = rtmp.org + rtmp.dir * ztmp;
        dv = glm::normalize(vtmp - v);
        ao += dot(n, dv); // divide by d at the end to normalize dot product
        nsamples++;
      }
      ao /= fmax(1, nsamples);
      imgrasters.colbuffer[y * w + x] =
          color_to_rgb(rgb_to_color(imgrasters.colbuffer[y * w + x]) *
                       (1 - glm::clamp<double>(ao, 0, 1)));
    }
  }
}

color ray_trace_pixel(double x, double y, int w, int h, const scene &sc,
                      const render_options &opts) {
  return trace_ray(sc.cam.cast_ray(x, y, w, h), sc, 1, opts.ray_depth);
}

color path_trace_pixel(double x, double y, int w, int h, const scene &sc,
                       const render_options &opts) {
  unsigned int s;
  color totalcol{};
  for (s = 1; s <= opts.samples; s++) {
    totalcol += trace_path(sc.cam.cast_ray(x, y, w, h), sc);
  }
  return tone_map(totalcol / static_cast<double>(s));
}

color amb_occ_pixel(double x, double y, int w, int h, const scene &sc,
                    const render_options &opts) {
  ray r = sc.cam.cast_ray(x, y, w, h);
  double ao = ambient_occlusion(r, sc, opts.samples);
  ao = glm::clamp<double>(ao * 2, 0, 1);
  return color(ao, ao, ao);
}

color ray_march_pixel(double x, double y, int w, int h, const scene &sc,
                      const render_options &opts) {
  int num_steps = 150;
  ray viewray = sc.cam.cast_ray(x, y, w, h);
  int steps;
  // TODO get this number
  //    double iterations;
  double t;
  if (!ray_march(viewray, sc.de_obj, &t, &steps, num_steps)) {
    return sc.w->get_color(viewray);
  } else {
    double ambient_occlusion = (1.0 - double(steps) / num_steps);
    //        double hue = 360 * iterations;
    return ambient_occlusion * color{1, 1, 1};
  }
}

__v4si edge_vect::init(const point_2d<int> &v0, const point_2d<int> &v1,
                       const point_2d<int> &origin) {
  // Edge setup
  int A = v0.y - v1.y, B = v1.x - v0.x;
  int C = v0.x * v1.y - v0.y * v1.x;

  // step deltas
  one_step_x = _mm_set1_epi32(A * step_x);
  one_step_y = _mm_set1_epi32(B * step_y);

  __v4si x = _mm_set1_epi32(origin.x) + _mm_set_epi32(3, 2, 1, 0);
  __v4si y = _mm_set1_epi32(origin.y);

  // barycentric coordinates at edges:
  __v4si out = muli32(_mm_set1_epi32(A), x) + muli32(_mm_set1_epi32(B), y) +
               _mm_set1_epi32(C);
  return out;
}
