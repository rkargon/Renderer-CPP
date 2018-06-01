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

  if (range_contains(mirror_illum, mat.illum)) {
    // TODO shininess
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

EmissionBSDF::EmissionBSDF(color c) : col(c){};

color EmissionBSDF::getLight(color, vertex, vertex, vertex) const {
  return col;
}

vertex EmissionBSDF::getIncidentDirection(vertex normal,
                                          vertex viewDirection) const {
  // does not actually matter, material does not reflect any incident light
  return glm::reflect(viewDirection, normal);
}

DiffuseBSDF::DiffuseBSDF(color c) : col(c){};

color DiffuseBSDF::getLight(color incidentColor, vertex incidentDirection,
                            vertex normal, vertex) const {
  return col * incidentColor * dot(normal, incidentDirection);
}

vertex DiffuseBSDF::getIncidentDirection(vertex normal,
                                         vertex viewDirection) const {
  vertex d;
  double cos = 1, r = 0;
  int i = 0;
  while (cos > r) {
    i++;
    d = random_direction();
    cos = glm::dot(normal, d);
    r = double(rand()) / double(RAND_MAX);
  }
  if (glm::dot(d, normal) * glm::dot(viewDirection, normal) > 0) {
    d *= -1;
  }
  return d;
}

GlossyBSDF::GlossyBSDF(color c, double r)
    : col(c), roughness(glm::clamp(r, 0.0, 1.0)){};

color GlossyBSDF::getLight(color incidentColor, vertex, vertex, vertex) const {
  return col * incidentColor;
}

vertex GlossyBSDF::getIncidentDirection(vertex normal,
                                        vertex viewDirection) const {
  // returns a random ray in a cone centered around the reflection of
  // viewDirection
  // the cone's width depends on the roughness
  return glm::normalize(glm::reflect(viewDirection, normal) +
                        random_direction() * roughness);
}

MixBSDF::MixBSDF(double f, BSDF *m1, BSDF *m2) : fac(f), mat1(m1), mat2(m2) {
  materialCounter = true;
};

color MixBSDF::getLight(color incidentColor, vertex incidentDirection,
                        vertex normal, vertex returningDirection) const {
  materialCounter = !materialCounter;
  return (materialCounter ? mat1 : mat2)
      ->getLight(incidentColor, incidentDirection, normal, returningDirection);
}

vertex MixBSDF::getIncidentDirection(vertex normal,
                                     vertex viewDirection) const {
  return (materialCounter ? mat1 : mat2)
      ->getIncidentDirection(normal, viewDirection);
}
