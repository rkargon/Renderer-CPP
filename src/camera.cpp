//
//  camera.cpp
//  Renderer
//
//  Created by Raphael Kargon on 6/8/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "camera.h"

camera::camera()
    : center(0), focus(0), normal(1, 0, 0), vert(0, 0, 1), fov(0.75),
      mindist(0.01), maxdist(100), ortho(false), dof_focus_distance(0),
      aperture_size(0) {
  calc_image_vectors();
}

/* Projections */

// Cast a ray from the screen into 3D space
ray camera::cast_ray(const double px, const double py, const double w,
                     const double h) const {
  double img_w = 2 * tan(fov / 2);
  double img_h = (img_w * h / w);
  double x = px / w - 0.5;
  double y = py / h - 0.5;
  ray castray;

  if (ortho) {
    castray.dir = normal;
    castray.org = center + _cx * (x * img_w) +
                  _cy * (-y * img_h); //-y because of screen coordinates
    return castray;
  } else {
    if (this->aperture_size == 0) {
      castray.org = center;
      castray.dir = normal + _cx * (x * img_w) + _cy * (-y * img_h);
    } else {
      auto dp = random_point_in_unit_circle();
      // TODO probably wrong :-(
      castray.org = center + (aperture_size * _cx * dp.x) +
                    (aperture_size * _cy * dp.y); // ray is shifted along image
                                                  // plane by random amount
                                                  // depending on aperture size
      vertex focal_point = center +
                           dof_focus_distance * (normal + _cx * (x * img_w) +
                                                 _cy * (-y * img_h));
      castray.dir = focal_point - castray.org;
    }
    normalize_in_place(castray.dir);
    return castray;
  }
}

// Project a vertex into 3D space onto a screen
point_2d<double> camera::project_vertex(const vertex &v, const double w,
                                        const double h) const {
  vertex dv = v - center;
  double x, y;

  // in perspective mode, only the relative angle (ie cosine, ie dot product)
  //  of the vector matters
  if (!ortho) {
    normalize_in_place(dv);
  }
  double dvdotnorm = glm::dot(dv, normal);
  if (dvdotnorm <= mindist || dvdotnorm > maxdist) {
    return point_2d<double>(nan(""), nan(""));
  }
  x = glm::dot(dv, _cx);
  y = glm::dot(dv, _cy);

  double px = (0.5 + x / fov) * w;
  double py = (0.5 - y * w / (h * fov)) * h;
  return point_2d<double>(px, py);
}

// The distance of a vertex from the camera.
double camera::vertex_depth(const vertex &v) const {
  vertex dv = v - center;
  if (ortho) {
    return glm::dot(dv, normal);
  } else {
    return glm::length(dv) * (glm::dot(dv, normal) < 0 ? -1 : 1);
  }
}

// Get vector from camera to given vertex
vertex camera::view_vector(const vertex &v) const {
  return ortho ? normal * glm::dot(v, normal) : v - center;
}

double camera::face_depth(const face &f) const {
  return vertex_depth(f.center());
}

/* manipulation of camera */

void camera::center_focus() {
  vertex camdist = focus - center;
  center = focus - normal * glm::length(camdist);
}

void camera::shift_focus(const double dx, const double dy) {
  focus += _cx * dx + _cy * dy;
}

void camera::zoom(const double zoomfactor) {
  center = focus + (center - focus) * zoomfactor;
}

vertex camera::image_plane_vector(const double dx, const double dy) {
  return _cx * dx + _cy * dy;
}

/* Rotation */
void camera::rotate_axis(const vertex &axis, const double dtheta) {
  vert = rotate(vert, axis, dtheta);
  normal = rotate(normal, axis, dtheta);
  calc_image_vectors();
}

void camera::rotate_local_x(const double dtheta) {
  rotate_axis(cross(normal, vert), dtheta);
}
void camera::rotate_local_y(const double dtheta) { rotate_axis(vert, dtheta); }
void camera::rotate_local_z(const double dtheta) {
  rotate_axis(normal, dtheta);
}

void camera::set_global_rotation(const double theta, const double rho,
                                 const double psi) {
  double st = std::sin(theta), ct = std::cos(theta), sr = std::sin(rho),
         cr = std::cos(rho), sp = std::sin(psi), cp = std::cos(psi);
  vert = vertex(-sr, cr * sp, cr * cp);
  normal = vertex(ct * cr, ct * sr * sp - cp * st, st * sp + ct * cp * sr);
  calc_image_vectors();
}

void camera::calc_image_vectors() {
  _cx = glm::normalize(glm::cross(normal, vert));
  _cy = glm::normalize(glm::cross(_cx, normal));
}
