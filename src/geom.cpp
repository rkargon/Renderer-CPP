//
//  geom.cpp
//  Renderer
//
//  Created by Raphael Kargon on 6/5/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "geom.h"

#include "mesh.h"

#include <glm/gtx/io.hpp>
#include <glm/gtx/norm.hpp>

/* VERTEX */

std::size_t std::hash<vertex>::operator()(const vertex &v) const noexcept {
  std::size_t seed = 0;
  for (int i = 0; i < 3; ++i) {
    hash_combine(seed, v[i]);
  }
  return seed;
}

/* FACE */
face::face() = default;

face::face(const vertex &norm, vertex_id v1, vertex_id v2, vertex_id v3,
           mesh *object, bool check_normal)
    : normal(norm), v{v1, v2, v3}, obj(object) {
  if (check_normal) {
    if (glm::length2(normal) > EPSILON && is_perpendicular(normal)) {
      normalize_in_place(normal);
    } else {
      normal = generate_normal(); // if normal is invalid, generate it
    }
  }
}

vertex &face::get_vert(unsigned int i) { return obj->vertices[v[i]]; }
const vertex &face::get_vert(unsigned int i) const {
  return obj->vertices[v[i]];
}

bool face::is_perpendicular(const vertex &normal) const {
  vertex v12 = get_vert(1) - get_vert(0);
  vertex v23 = get_vert(2) - get_vert(1);
  // normal must be perpendicular to two of the edges, and must have nonzero
  // length
  return (!eps_zero(glm::length2(normal)) && eps_zero(glm::dot(v12, normal)) &&
          eps_zero(glm::dot(v23, normal)));
}
// use cross product of edges to create normal. Orientation depends on order of
// vertices
vertex face::generate_normal() const {
  return face::generate_normal(get_vert(0), get_vert(1), get_vert(2));
}

vertex face::generate_normal(const vertex &v1, const vertex &v2,
                             const vertex &v3) {
  vertex v12 = v2 - v1;
  vertex v23 = v3 - v2;
  vertex n = glm::cross(v12, v23);
  return glm::normalize(n);
}

// the center of the face, the arithmetic mean of vertices.
vertex face::center() const {
  return (get_vert(0) + get_vert(1) + get_vert(2)) * (1 / 3.0);
}

// tuv stores:
// t: distance from origin to intersect point
// u, v: barycentric coordinates of intersect based on edge1, edge2
// Muller-Trumbore algorithm
// TODO return optional?
bool face::intersect_ray_triangle(const ray &r, vertex *tuv) const {
  double t, u, v;

  vertex edge1 = get_vert(1) - get_vert(0);
  vertex edge2 = get_vert(2) - get_vert(0);

  vertex pvec = glm::cross(r.dir, edge2);
  double det = glm::dot(edge1, pvec);

  vertex tvec = r.org - get_vert(0);
  double inv_det = 1.0 / det;
  vertex qvec = glm::cross(tvec, edge1);

  if (det > EPSILON) {
    u = glm::dot(tvec, pvec);
    if (u < 0.0 || u > det) {
      return false;
    }
    v = glm::dot(r.dir, qvec);
    if (v < 0.0 || u + v > det) {
      return false;
    }
  } else if (det < -EPSILON) {
    u = glm::dot(tvec, pvec);
    if (u > 0.0 || u < det) {
      return false;
    }
    v = glm::dot(r.dir, qvec);
    if (v > 0.0 || u + v < det) {
      return false;
    }
  } else {
    return false;
  }

  u *= inv_det;
  v *= inv_det;
  t = glm::dot(edge2, qvec) * inv_det;
  if (tuv != nullptr) {
    *tuv = vertex(t, u, v);
  }
  return (t > EPSILON);
}

bounds face::bounding_box() const {
  bounds b;
  b.min = vertex(min_coord(0), min_coord(1), min_coord(2));
  b.max = vertex(max_coord(0), max_coord(1), max_coord(2));
  return b;
}

std::ostream &operator<<(std::ostream &os, const face &f) {
  os << "face: n=" << f.normal << ", obj=" << f.obj << ", (" << f.v[0] << ", "
     << f.v[1] << ", " << f.v[2] << ")" << std::endl;
  for (int i = 0; i < 3; ++i) {
    os << "f.v" << i << ": " << f.get_vert(i) << std::endl;
  }
  return os;
}

bounds::bounds() : min(0, 0, 0), max(0, 0, 0) {}
bounds::bounds(const vertex &a, const vertex &b) : min(a), max(b) {}

// surface area
double bounds::area() const {
  double dx = fabs(max.x - min.x);
  double dy = fabs(max.y - min.y);
  double dz = fabs(max.z - min.z);
  return 2 * (dx * dy + dy * dz + dz * dx);
}

std::ostream &operator<<(std::ostream &os, bounds &b) {
  return os << b.min << " --> " << b.max;
}

/*  EDGE  */

edge::edge(vertex_id v1, vertex_id v2) { this->set(v1, v2); }

bool edge::operator==(const edge &e2) const {
  return (this->v1 == e2.v1) && (this->v2 == e2.v2);
}

void edge::set(vertex_id v1, vertex_id v2) {
  std::tie(this->v1, this->v2) = std::minmax(v1, v2);
}

const vertex_id &edge::least() const { return v1; }
const vertex_id &edge::greatest() const { return v2; }
const vertex_id &edge::first() const { return v1; }
const vertex_id &edge::second() const { return v2; }
const vertex_id &edge::operator[](std::size_t i) const { return this->vec[i]; }

std::size_t std::hash<edge>::operator()(const edge &e) const noexcept {
  std::size_t hash_1 = std::hash<vertex_id>{}(e.least());
  hash_combine(hash_1, e.greatest());
  return hash_1;
}

std::ostream &operator<<(std::ostream &os, const edge &e) {
  return os << "edge(" << e.first() << ", " << e.second() << ")";
}

/*  MISC  */

// Rotates the given vertex about an axis, the given number of degrees
vertex rotate(const vertex &v, const vertex &axis, const double dtheta) {
  vertex a = glm::normalize(axis);
  double l = a.x, m = a.y, n = a.z;
  double s = std::sin(dtheta), c = std::cos(dtheta);

  // columns of rotation matrix
  vertex col1(l * l * (1 - c) + c, l * m * (1 - c) + n * s,
              l * n * (1 - c) - m * s);
  vertex col2(m * l * (1 - c) - n * s, m * m * (1 - c) + c,
              m * n * (1 - c) + l * s);
  vertex col3(n * l * (1 - c) + m * s, n * m * (1 - c) - l * s,
              n * n * (1 - c) + c);

  return vertex(glm::dot(v, col1), glm::dot(v, col2), glm::dot(v, col3));
}

// A random unit vector.
vertex random_direction() {
  double theta = (double)rand() / RAND_MAX;
  theta *= M_PI * 2;
  double z = (double)rand() / RAND_MAX;
  z = z * 2 - 1;
  double w = std::sqrt(1 - z * z);
  return vertex(w * std::cos(theta), w * std::sin(theta), z);
}

vertex lerp(const vertex &v1, const vertex &v2, const vertex &v3, double w1,
            double w2, double w3) {
  return v1 * w1 + v2 * w2 + v3 * w3;
}
vertex lerp(const vertex &v1, const vertex &v2, const vertex &v3,
            const vertex &w) {
  return v1 * w.x + v2 * w.y + v3 * w.z;
}

vertex from_polar(const double radius, const double theta, const double phi) {
  return radius * vertex(std::cos(phi) * std::sin(theta),
                         std::sin(phi) * std::sin(theta), std::cos(theta));
}

void to_polar(const vertex &v, double &radius, double &theta, double &phi) {
  radius = glm::length(v);
  theta = std::acos(v.z / radius);
  phi = std::atan2(v.y, v.x);
}

void to_polar_angles(const vertex &v, const double radius, double &theta,
                     double &phi) {
  theta = std::acos(v.z / radius);
  phi = std::atan2(v.y, v.x);
}

vertex box_fold(const vertex &v, const double l) {
  return 2.0 * glm::clamp(v, vertex(-l), vertex(l)) - v;
}

// TODO test if optimization makes a difference
double sphere_fold_ratio(const vertex &v, const double min_radius_2,
                         const double fixed_radius_2) {
  double r2 = glm::length2(v);
  if (r2 < min_radius_2) {
    return fixed_radius_2 / min_radius_2;
  } else if (r2 < fixed_radius_2) {
    return fixed_radius_2 / r2;
  } else {
    return 1;
  }

  //    double rad = this->len();
  //    /*  branchless optimization */
  //    double f = (rad < r) * ((r * r) / (rad * rad)) + (rad >= r) *
  //     ((rad < 1) / rad + (rad >= 1)); return *this * f;
}

// color functions

double luma(const color &c) {
  return glm::dot(c, color(0.2126, 0.7152, 0.0722));
}

// Take a vertex of three values in range [0,1] and convert to 32-bit RGB
uint color_to_rgb(const color &c) {
  return (unsigned char)glm::clamp(c.r * 255.0, 0.0, 255.0) << 16 |
         (unsigned char)glm::clamp(c.g * 255.0, 0.0, 255.0) << 8 |
         (unsigned char)glm::clamp(c.b * 255.0, 0.0, 255.0);
}
// Convert a normal vector to 32-bit RGB. FOr each coordinate, the range [-1, 1]
// is mapped to [0.0, 255]
uint normal_to_rgb(const vertex &n) {
  return (unsigned char)glm::clamp(n.x * 128 + 128, 0.0, 255.0) << 16 |
         (unsigned char)glm::clamp(n.y * 128 + 128, 0.0, 255.0) << 8 |
         (unsigned char)glm::clamp(n.z * 128 + 128, 0.0, 255.0);
}
// Convert 32-bit RGB to vertex
color rgb_to_color(const uint rgb) {
  return color(((rgb >> 16) & 0xff) / 255.0, ((rgb >> 8) & 0xff) / 255.0,
               (rgb & 0xff) / 255.0);
}
// Convert 32-bit RGB to normal vector
color rgb_to_normal(const uint n) {
  return vertex(((n >> 16) & 0xff) / 255.0 - 0.5,
                ((n >> 8) & 0xff) / 255.0 - 0.5, (n & 0xff) / 255.0 - 0.5);
}

color hsv_to_rgb(const int hue, const double saturation, const double value) {
  double chroma = value * saturation;
  double hue_mod = hue / 60.0;
  double x = chroma * (1 - fabs(fmod(hue_mod, 2) - 1));
  switch (int(hue_mod)) {
  case 0:
    return {chroma, x, 0};
  case 1:
    return {x, chroma, 0};
  case 2:
    return {0, chroma, x};
  case 3:
    return {0, x, chroma};
  case 4:
    return {x, 0, chroma};
  case 5:
    return {chroma, 0, x};
  default:
    return {0, 0, 0};
  }
}

color tone_map(const color &c) {
  double l = luma(c);
  return c / (l + 1);
}

/* bounding box related functions */

// THe bounding box of a set of faces.
bounds calc_bounding_box(const std::vector<const face *> &faces) {
  double minx, miny, minz, maxx, maxy, maxz, tmp;
  minx = miny = minz = maxx = maxy = maxz = nan("");
  for (const face *f : faces) {
    tmp = f->min_coord(0);
    if (isnan(minx) || minx > tmp) {
      minx = tmp;
    }
    tmp = f->min_coord(1);
    if (isnan(miny) || miny > tmp) {
      miny = tmp;
    }
    tmp = f->min_coord(2);
    if (isnan(minz) || minz > tmp) {
      minz = tmp;
    }

    tmp = f->max_coord(0);
    if (isnan(maxx) || maxx < tmp) {
      maxx = tmp;
    }
    tmp = f->max_coord(1);
    if (isnan(maxy) || maxy < tmp) {
      maxy = tmp;
    }
    tmp = f->max_coord(2);
    if (isnan(maxz) || maxz < tmp) {
      maxz = tmp;
    }
  }
  bounds b;
  b.min = vertex(minx, miny, minz);
  b.max = vertex(maxx, maxy, maxz);
  return b;
}
// Stores intersection of b1 and b2 in newbounds
void intersect_bounding_boxes(const bounds &b1, const bounds &b2,
                              bounds &newbounds) {
  newbounds.min = max(b1.min, b2.min);
  newbounds.max = min(b1.max, b2.max);
}
void list_bounds(bounds &newbounds, int nverts, const vertex *v1...) {
  newbounds.min = vertex(NAN, NAN, NAN);
  newbounds.max = vertex(NAN, NAN, NAN);
  const vertex *vtmp;

  va_list vertices;
  va_start(vertices, v1);
  vtmp = v1;
  for (int i = 0; i < nverts; i++) {
    if (newbounds.min != newbounds.min || newbounds.min.x > vtmp->x)
      newbounds.min.x = vtmp->x;
    if (newbounds.min != newbounds.min || newbounds.min.y > vtmp->y)
      newbounds.min.y = vtmp->y;
    if (newbounds.min != newbounds.min || newbounds.min.z > vtmp->z)
      newbounds.min.z = vtmp->z;
    if (newbounds.max != newbounds.max || newbounds.max.x < vtmp->x)
      newbounds.max.x = vtmp->x;
    if (newbounds.max != newbounds.max || newbounds.max.y < vtmp->y)
      newbounds.max.y = vtmp->y;
    if (newbounds.max != newbounds.max || newbounds.max.z < vtmp->z)
      newbounds.max.z = vtmp->z;
    vtmp = va_arg(vertices, vertex *);
  }
}

// geometry intersection

// Intersect a ray with a bounding box (AABB = Axis Aligned Bounding Box)
bool ray_AABB_intersect(const bounds &AABB, const ray &r) {
  double tmp;
  double tmin = (AABB.min.x - r.org.x) / r.dir.x;
  double tmax = (AABB.max.x - r.org.x) / r.dir.x;
  if (tmin > tmax) {
    tmp = tmin;
    tmin = tmax;
    tmax = tmp;
  }

  double tymin = (AABB.min.y - r.org.y) / r.dir.y;
  double tymax = (AABB.max.y - r.org.y) / r.dir.y;
  if (tymin > tymax) {
    tmp = tymin;
    tymin = tymax;
    tymax = tmp;
  }
  if (tmin > tymax || tmax < tymin)
    return false;
  if (tymin > tmin)
    tmin = tymin;
  if (tymax < tmax)
    tmax = tymax;

  double tzmin = (AABB.min.z - r.org.z) / r.dir.z;
  double tzmax = (AABB.max.z - r.org.z) / r.dir.z;
  if (tzmin > tzmax) {
    tmp = tzmin;
    tzmin = tzmax;
    tzmax = tmp;
  }
  if (tmin > tzmax || tmax < tzmin)
    return false;
  if (tzmin > tmin)
    tmin = tzmin;
  if (tzmax < tmax)
    tmax = tzmax;

  if (tmin <= 0 && tmax <= 0)
    return false; // ray should not intersect bounding box if box is behind the
                  // origin of the ray
  return true;
}
// Intersect a ray with a set of faces. If lazy==false, return nearest face.
// Otherwise, return first intersection.
//(t,u,v) stores (intersection distance, barycentric coordinate from v12,
// barycentric coordiante from v13)
const face *ray_faces_intersect(const std::vector<const face *> &faces,
                                const ray &r, bool lazy, vertex *tuv) {
  const face *f = nullptr;
  vertex tuvtmp;
  double zmin = std::nan("");
  for (const face *ftmp : faces) {
    if (ftmp->intersect_ray_triangle(r, &tuvtmp) &&
        (tuvtmp[0] < zmin || std::isnan(zmin))) {
      zmin = tuvtmp[0];
      f = ftmp;
      if (tuv != nullptr) {
        *tuv = tuvtmp;
      }
      if (lazy) {
        return f;
      }
    }
  }
  return f;
}
// Intersect a ray with a sphere, t stores distance from ray origin to intersect
bool ray_sphere_intersect(const ray &r, const double rad, double &t) {
  double d_dot_o = dot(r.dir, r.org);
  double lensq = glm::length2(r.org);
  double disc = d_dot_o * d_dot_o - (lensq - rad * rad);
  if (disc < 0) {
    return false;         // no intersect
  } else if (disc == 0) { // one intersect
    t = -d_dot_o;
    return (t >= 0);
  } else {
    // return nearest intersection, that is not behind ray origin
    disc = std::sqrt(disc);
    // two intersection points.
    double t1 = disc - d_dot_o, t2 = (-disc) - d_dot_o;
    if (t1 >= 0 && t2 >= 0) {
      t = std::fmin(t1, t2);
      return true;
    } else if (t1 >= 0 || t2 >= 0) {
      t = std::fmax(t1, t2);
      return true;
    } else
      return false;
  }
}

point_2d<double> random_point_in_unit_circle() {
  double theta = double(rand()) / RAND_MAX * M_PI * 2;
  double r = double(rand()) / RAND_MAX;
  return {r * cos(theta), r * sin(theta)};
}
