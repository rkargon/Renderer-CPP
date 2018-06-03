//
//  BSDF.cpp
//  Renderer
//
//  Created by Raphael Kargon on 11/9/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "BSDF.h"

#include "glm/gtc/type_ptr.hpp"
#include "glm/gtx/component_wise.hpp"

#include "tiny_obj_loader/tiny_obj_loader.h"

std::unique_ptr<BSDF> BSDF::from_obj_mat(const tinyobj::material_t &mat) {
  int mirror_illum[] = {3, 5, 8};
  int glass_illum[] = {4, 6, 7, 9};

  /// TODO TEMPORARY XXX

  if (range_contains(mirror_illum, mat.illum)) {
    // TODO shininess
    return std::make_unique<DiffuseBSDF>(glm::make_vec3(mat.diffuse));
    return std::make_unique<GlossyBSDF>(
        glm::make_vec3(mat.specular), 0 * 1.0 / std::fmax(1.0, mat.shininess));
  } else if (range_contains(glass_illum, mat.illum)) {
    // TODO GLASS
    return std::make_unique<DiffuseBSDF>();
  } else {
    color emission = glm::make_vec3(mat.emission);
    if (glm::compMax(emission) > 0) {
      return std::make_unique<EmissionBSDF>(emission);
    } else {
      return std::make_unique<DiffuseBSDF>(glm::make_vec3(mat.diffuse));
    }
  }
}

DiffuseBSDF::DiffuseBSDF(color c) : col(c){};

direction_sample
DiffuseBSDF::sample_direction(const vertex &normal,
                              const vertex &outgoing_direction) const {
  // TODO use importance sampling
  return uniform_unit_hemisphere(normal);
}

color DiffuseBSDF::bsdf(const vertex &incoming_dir,
                        const vertex &outgoing_dir) const {
  return col / M_PI;
}

EmissionBSDF::EmissionBSDF(color c) : col(c){};

direction_sample
EmissionBSDF::sample_direction(const vertex &normal,
                               const vertex &outgoing_direction) const {
  // doesn't really matter, emitters don't reflect incoming light.
  return {normal, 1.0};
}

color EmissionBSDF::bsdf(const vertex &incoming_dir,
                         const vertex &outgoing_dir) const {
  return color(1.0);
}

GlossyBSDF::GlossyBSDF(color c, double r)
    : col(c), roughness(glm::clamp(r, 0.0, 1.0)){};

direction_sample
GlossyBSDF::sample_direction(const vertex &normal,
                             const vertex &outgoing_direction) const {
  // TODO importance sampling
  return uniform_unit_hemisphere(normal);
}

color GlossyBSDF::bsdf(const vertex &incoming_dir,
                       const vertex &outgoing_dir) const {
  // TODO
  return color(1.0);
}

MixBSDF::MixBSDF(double f, BSDF *m1, BSDF *m2) : fac(f), mat1(m1), mat2(m2) {
  materialCounter = true;
};

direction_sample
MixBSDF::sample_direction(const vertex &normal,
                          const vertex &outgoing_direction) const {
  // TODO importance sampling
  return uniform_unit_hemisphere(normal);
}

color MixBSDF::bsdf(const vertex &incoming_dir,
                    const vertex &outgoing_dir) const {
  // TODO
  return color(1.0);
}
