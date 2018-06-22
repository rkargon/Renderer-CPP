//
//  scene.cpp
//  Renderer
//
//  Created by Raphael Kargon on 6/8/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "scene.h"

#include "mesh.h"

#include "glm/gtc/type_ptr.hpp"
#include "glm/gtx/component_wise.hpp"
#include "glm/gtx/io.hpp"
#include "tiny_obj_loader/tiny_obj_loader.h"

#include <iostream>
#include <unordered_map>
#include <unordered_set>

material::material()
    : diff_col(0), spec_col(0), spec_hardness(128), refl_intensity(0), alpha(1),
      ior(1), col_tex(nullptr), spec_tex(nullptr), norm_tex(nullptr) {}

material::material(const color &diff_col) : material() {
  this->diff_col = diff_col;
}

material::material(const tinyobj::material_t &mat)
    : diff_col(glm::make_vec3(mat.diffuse)),
      spec_col(glm::make_vec3(mat.specular)), spec_hardness(mat.shininess),
      refl_intensity(static_cast<double>(mat.is_reflective())),
      alpha(1.0 - glm::compMax(glm::make_vec3(mat.transmittance))),
      ior(mat.ior), col_tex(nullptr), spec_tex(nullptr), norm_tex(nullptr) {}

// Either returns diff_col, or if col_tex is defined, the pixel corresponding to
// the passed texture coordinates
color material::get_color(double txu, double txv) const {
  if (col_tex == nullptr)
    return diff_col;
  else
    return rgb_to_color(
        col_tex->pixel(txu * col_tex->width(), txv * col_tex->height()));
}
// TODO: Return normal relative to a passed, default normal
vertex material::get_normal(double txu, double txv) const {
  if (norm_tex == nullptr)
    return vertex();
  else
    return rgb_to_color(
        norm_tex->pixel(txu * norm_tex->width(), txv * norm_tex->height()));
}
// Returns spec_col, or the corresponding pixel of spec_tex to determine
// specular color
color material::get_spec_col(double txu, double txv) const {
  if (spec_tex == nullptr)
    return spec_col;
  else
    return rgb_to_color(
        spec_tex->pixel(txu * spec_tex->width(), txv * spec_tex->height()));
}

lamp::lamp() : intensity(1), falloff(2), loc(vertex()), col(color(1, 1, 1)) {}

lamp::lamp(const double i, const int fall, const vertex &l, const color &c)
    : intensity(i), falloff(fall), loc(l), col(c), is_sun(false) {}

lamp::lamp(const double i, const int fall, const vertex &l, const color &c,
           const vertex &sun_direction)
    : intensity(i), falloff(fall), loc(l), col(c), is_sun(true),
      sun_direction(sun_direction) {}

world::world()
    : horizon_col(0), zenith_col(0), is_flat(true), ambient_intensity(0) {}

// PRECONDITION: dir is normalized
// If is_flat, horizon_col is returned.
// Otherwise, z component of r.z is used to interpolate between horizon_col and
// zenith_col
color world::get_color(const ray &r) const {
  if (is_flat)
    return horizon_col;
  else {
    return glm::mix(horizon_col, zenith_col,
                    fabs(r.dir.z)); // not really a linear progression
                                    // from horizon to zenith, but
                                    // close enough
  }
}

// sky::sky(vertex beta_r, vertex beta_m, double Hr, double Hm,
//          double radius_earth, double radius_atmo, vertex sun_direction,
//          double sun_intensity, double g)
//     : beta_r(beta_r), beta_m(beta_m), Hr(Hr), Hm(Hm),
//       radius_earth(radius_earth), radius_atmo(radius_atmo),
//       sun_direction(glm::normalize(sun_direction)),
//       sun_intensity(sun_intensity), g(g) {}

sky::sky()
    : beta_r(5.5e-6, 13.0e-6, 22.4e-6), beta_m(21e-6, 21e-6, 21e-6), Hr(7994),
      Hm(1200), radius_earth(6360e3), radius_atmo(6420e3),
      sun_direction(glm::normalize(vertex(0, -1, 0.4))), sun_intensity(20),
      g(0.76) {}

// Source:
// http://scratchapixel.com/lessons/3d-advanced-lessons/simulating-the-colors-of-the-sky/atmospheric-scattering/
// Given aa view ray, calculates the color of sky in that direction using
// atmospheric scattering.
color sky::get_color(const ray &r) const {
  ray r_new = r;
  r_new.org.z += radius_earth;
  // get atmosphere intersection point
  double t_atmo;
  if (!ray_sphere_intersect(r_new, radius_atmo, t_atmo)) {
    return color(0, 0, 0);
  }
  double segment_length = t_atmo / samples;
  vertex ray_segment = r_new.dir * segment_length,
         ray_sample = r_new.org + ray_segment * 0.5; // set up initial sample
                                                     // and delta for each
                                                     // sample along ray

  // get rayleigh and mie phase functions
  double mu = dot(r_new.dir, sun_direction);
  double phase_r = 3 * (1 + mu * mu) / (16 * M_PI);
  double phase_m = 3 * (1 - g * g) * (1 + mu * mu) /
                   (8 * M_PI * (2 + g * g) * pow(1 + g * g - 2 * g * mu, 1.5));

  double optical_depth_r = 0, optical_depth_m = 0;
  color col_r(0), col_m(0); // resulting colors of each type of scattering
  // for each sample along ray
  for (int i = 0; i < samples; i++, ray_sample += ray_segment) {
    // get sample altitude
    double height = glm::length(ray_sample) - radius_earth;
    // update optical depth
    double hr = std::exp(-height / Hr) * segment_length;
    double hm = std::exp(-height / Hm) * segment_length;
    optical_depth_r += hr;
    optical_depth_m += hm;

    // calc amount of sunlight coming along ray at sample point
    ray lightray(ray_sample, sun_direction);
    double lightray_t;
    // cast ray from sample to sun
    ray_sphere_intersect(lightray, radius_atmo, lightray_t);
    double segment_length_light = lightray_t / samples_lightray;
    vertex lightray_segment = lightray.dir * segment_length_light;
    vertex lightray_sample = lightray.org + lightray_segment * 0.5;
    double optical_depth_light_r = 0, optical_depth_light_m = 0;

    // test several points along this ray
    int j;
    for (j = 0; j < samples_lightray;
         j++, lightray_sample += lightray_segment) {
      double heightlight = glm::length(lightray_sample) - radius_earth;
      if (heightlight < 0) {
        break; // discard ray if it is in shadow of the earth
      }
      // calc optical depth at each subsample
      optical_depth_light_r +=
          std::exp(-heightlight / Hr) * segment_length_light;
      optical_depth_light_m +=
          std::exp(-heightlight / Hm) * segment_length_light;
    }
    if (j == samples_lightray) { // ray is not in shadow of earth
      // calc light color based on attenuation
      vertex tau = beta_r * (optical_depth_r + optical_depth_light_r) +
                   beta_m * 1.1 * (optical_depth_m + optical_depth_light_m);
      vertex attenuation(std::exp(-tau.x), std::exp(-tau.y), std::exp(-tau.z));
      col_r += hr * attenuation;
      col_m += hm * attenuation;
    } else {
      return color();
    }
  }
  return (col_r * beta_r * phase_r + col_m * beta_m * phase_m) * sun_intensity;
}

scene::scene() : w(std::make_unique<world>()) {}

scene::scene(const std::string &filename) : scene() {
  tinyobj::attrib_t attrs;
  std::vector<tinyobj::shape_t> shapes;
  std::vector<tinyobj::material_t> mats;
  std::string err;
  // TODO use stdx::filesystem
  // TODO check if empty material basedir works
  PRINT_VARIABLE(filename);
  PRINT_VARIABLE(containing_directory(filename));
  bool succ = tinyobj::LoadObj(&attrs, &shapes, &mats, &err, filename.c_str(),
                               containing_directory(filename).c_str(), true);
  if (!succ) {
    throw std::runtime_error("Failed to load .obj file: " + err);
  }

  // convert OBJ materials
  if (mats.size() > 0) {
    this->materials.insert(this->materials.begin(), mats.begin(), mats.end());
    this->bsdfs.reserve(materials.size());
    std::transform(mats.begin(), mats.end(), std::back_inserter(this->bsdfs),
                   BSDF::from_obj_mat);
  } else {
    this->materials.emplace_back();
    this->bsdfs.push_back(std::make_unique<DiffuseBSDF>());
  }
  for (const auto &s : shapes) {
    objects.emplace_back(std::make_unique<mesh>());
    mesh &m = *objects.back();
    m.name = s.name;
    if (mats.size() > 0) {
      // TODO only one material per object right now.
      auto mat_id = s.mesh.material_ids[0];
      m.mat_id = mat_id;
      m.bsdf = this->bsdfs[mat_id].get();
      std::cout << "!!! " << mats[mat_id].name << ", mesh name: " << m.name
                << " diffuse: " << glm::make_vec3(mats[mat_id].diffuse)
                << std::endl;
    } else {
      m.mat_id = 0;
      m.bsdf = this->bsdfs.front().get();
    }

    m.faces.reserve(s.mesh.indices.size() / 3);
    // TODO this should be a vector
    std::unordered_map<vertex_id, vertex_id> vertex_ids_map; // old to new
    std::unordered_set<edge> edges;
    for (std::size_t i = 0; i < s.mesh.indices.size(); i += 3) {
      // get vertices & their IDs
      vertex v[3];
      vertex_id vid[3]; // new ids
      for (int vi = 0; vi < 3; ++vi) {
        auto iter_bool = vertex_ids_map.emplace(
            s.mesh.indices[i + vi].vertex_index, vertex_ids_map.size());
        vid[vi] = iter_bool.first->second;
        v[vi] = glm::make_vec3(
            &attrs.vertices[3 * s.mesh.indices[i + vi].vertex_index]);
      }
      // face
      vertex norm = face::generate_normal(v[0], v[1], v[2]);
      m.faces.push_back(face(norm, vid[0], vid[1], vid[2], &m, false));

      // add edges
      edges.emplace(vid[0], vid[1]);
      edges.emplace(vid[1], vid[2]);
      edges.emplace(vid[2], vid[0]);
    }

    // add vertices from hash
    m.vertices.resize(vertex_ids_map.size());
    for (const auto &kv : vertex_ids_map) {
      m.vertices[kv.second] = glm::make_vec3(&attrs.vertices[3 * kv.first]);
    }

    // add edges from hash
    m.edges.insert(m.edges.begin(), edges.begin(), edges.end());

    // face adjacencies
    m.face_adjacencies.resize(m.vertices.size());
    for (face_id fi = 0; fi < m.faces.size(); ++fi) {
      for (vertex_id vi : m.faces[fi].v) {
        m.face_adjacencies[vi].insert(fi);
      }
    }

    // calc vertex normals
    m.vertex_normals.resize(m.vertices.size());
    for (const auto &kv : vertex_ids_map) {
      vertex_id old_id = kv.first;
      vertex_id new_id = kv.second;
      vertex normal;
      if (attrs.normals.size() > 0) {
        normal = glm::make_vec3(&attrs.normals[3 * old_id]);
      } else {
        normal = vertex(0);
      }
      if (glm::length2(normal) > EPSILON) {
        normalize_in_place(normal);
      } else {
        vertex n(0);
        for (face_id f_id : m.face_adjacencies[new_id]) {
          n += m.faces[f_id].normal;
        }
        normal = glm::normalize(n);
      }

      m.vertex_normals[new_id] = normal;
    }
  }
}

mesh &scene::add_object(std::ifstream &infile, const std::string &name,
                        bool update_tree) {
  objects.push_back(std::make_unique<mesh>(infile, name));
  if (update_tree) {
    this->update_tree();
  }
  return *objects.back();
}

void scene::update_tree() {
  std::vector<const face *> allfaces;
  for (const auto &o_ptr : objects) {
    allfaces.reserve(allfaces.size() + o_ptr->faces.size());
    for (const face &f : o_ptr->faces) {
      allfaces.push_back(&f);
    }
  }
  this->kdt = kdtree::build_tree(allfaces);
}
