//
//  geom.h
//  Renderer
//
// Functions and data structures relating to geometry.
//
//  Created by Raphael Kargon on 6/5/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__geom__
#define __Renderer__geom__

#include "common.h"

#include "glm/glm.hpp"

#include <iostream>
#include <vector>

typedef glm::dvec3 vertex;
typedef vertex color;
typedef unsigned int vertex_id;
typedef unsigned int face_id;

class face;
class edge;
class mesh;

template <typename T> struct point_2d {
  T x, y;
  point_2d() {}
  point_2d(T a, T b) : x(a), y(b) {}
};
// cross product of AB, and AC
template <typename T> T orient_2d(point_2d<T> a, point_2d<T> b, point_2d<T> c) {
  return ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
}

// // Represents a point in three-dimensional space in cartesian coordinates.
// Also
// // used to represent color values.
// class vertex {
// public:
//   vertex();
//   vertex(double a, double b, double c);

//   inline double len() const { return sqrt(x * x + y * y + z * z); }
//   inline double lensqr() const { return (x * x + y * y + z * z); }
//   void normalize();
//   void set(const double a, const double b, const double c);
//   vertex reflection(const vertex &n) const;
//   vertex unitvect() const; // return a normalized copy of this vector
//   void to_polar(double &radius, double &theta, double &phi) const;
//   // Convenience function for getting polar coordinates, given a
//   pre-calculated
//   // radius. This saves some computation time.
//   void to_polar_angles(const double radius_squared, double &theta,
//                        double &phi) const;
//   vertex box_fold(const double l) const;
//   // Returns the scaling factor for a sphere fold
//   double sphere_fold_ratio(const double min_radius_2,
//                            const double fixed_radius_2) const;

//   void operator+=(const vertex &v2);
//   void operator+=(const vertex *v2);
//   void operator-=(const vertex &v2);
//   void operator-=(const vertex *v2);
//   void operator*=(const double r);
//   void operator*=(const vertex &v2);
//   void operator*=(const vertex *v2);
//   vertex operator-();
//   friend vertex operator+(const vertex &v1, const vertex &v2);
//   friend vertex operator+(const vertex &v1, const vertex *v2);
//   friend vertex operator-(const vertex &v1, const vertex &v2);
//   friend vertex operator-(const vertex &v1, const vertex *v2);
//   friend vertex operator*(const vertex &v1, double r);
//   friend vertex operator*(const double r, const vertex &v1);
//   friend vertex operator*(const vertex &v1, const vertex &v2);
//   friend vertex operator/(const double r, const vertex &v1);
//   friend vertex operator/(const vertex &v1, const double r);
//   friend bool operator==(const vertex &v1, const vertex &v2);
//   friend bool operator!=(const vertex &v1, const vertex &v2);
//   friend std::ostream &operator<<(std::ostream &os, const vertex &v);

//   // data
//   union {
//     struct {
//       double r, g, b;
//     };
//     struct {
//       double t, u, v;
//     };
//     struct {
//       double x, y, z;
//     };
//     struct {
//       double vec[3];
//     };
//     // ok so maybe I'm lazy
//   };
// };

template <typename GLM_T> void normalize_in_place(GLM_T &x) {
  x = glm::normalize(x);
}

// Represents an axis-aligned bounding box. Also represents rays, with an origin
// and direction
class bounds {
public:
  union {
    struct {
      vertex min, max;
    };
    struct {
      vertex org, dir;
    };
    struct {
      vertex varr[2];
    };
  };

  bounds();
  bounds(const vertex &a, const vertex &b);

  double area() const;
  inline double volume() const {
    return (max.x - min.x) * (max.y - min.y) * (max.z - min.z);
  }
  inline double d(const int k) const {
    return max[k] - min[k];
  } // doesn't check axis indices for validity. Be careful.
};

typedef bounds ray;

// Represents a triangular face on a 3D mesh
class face {
public:
  vertex normal;
  vertex_id v[3];
  mesh *obj;

  face();
  face(const vertex &norm, vertex_id v1, vertex_id v2, vertex_id v3,
       mesh *object);

  vertex &get_vert(unsigned int i);
  const vertex &get_vert(unsigned int i) const;

  bounds bounding_box() const;
  vertex center() const;
  vertex generate_normal();
  bool intersect_ray_triangle(const ray &r, vertex *tuv);
  bool is_perpendicular(const vertex &normal);
  inline double min_coord(int axis) const {
    return fmin(fmin(get_vert(0)[axis], get_vert(1)[axis]), get_vert(2)[axis]);
  }
  inline double max_coord(int axis) const {
    return fmax(fmax(get_vert(0)[axis], get_vert(1)[axis]), get_vert(2)[axis]);
  }
};

namespace std {
template <> struct hash<vertex> {
  std::size_t operator()(const vertex &v) const noexcept;
};
}
// an undirected edge. The hash for this ignores the ordering of the vertices.
// It is "undirected" for equality & hashing purposes, because during
// construction thw two vertices are sorted least to greatest
class edge {
public:
  edge(vertex_id v1, vertex_id v2);
  bool operator==(const edge &e2) const;
  void set(vertex_id v1, vertex_id v2);

  // getters, w/ various names for convenience
  const vertex_id &least() const;
  const vertex_id &greatest() const;
  const vertex_id &first() const;
  const vertex_id &second() const;

private:
  union {
    struct {
      vertex_id v1, v2;
    };
    vertex_id vec[2];
  };
};

namespace std {
template <> struct hash<edge> {
  std::size_t operator()(const edge &e) const noexcept;
};
}

std::ostream &operator<<(std::ostream &os, const edge &e);

/* vertex functions */
vertex rotate(const vertex &v, const vertex &axis, const double dtheta);
vertex random_direction(); // picks a random point on the unit sphere
                           // (uniformly)

// linearly interpolate between 3 vertices using barycentric coordinates
vertex lerp(const vertex &v1, const vertex &v2, const vertex &v3, double w1,
            double w2, double w3);
vertex lerp(const vertex &v1, const vertex &v2, const vertex &v3,
            const vertex &w);

vertex from_polar(const double radius, const double theta, const double phi);
void to_polar(const vertex &v, double &radius, double &theta, double &phi);
// Convenience function for getting polar coordinates, given a pre -
// calculated radius. This saves some computation time.
void to_polar_angles(const vertex &v, const double radius_squared,
                     double &theta, double &phi);

vertex box_fold(const vertex &v, const double l);
// Returns the scaling factor for a sphere fold
double sphere_fold_ratio(const vertex &v, const double min_radius_2,
                         const double fixed_radius_2);

/* color functions */
uint color_to_rgb(const color &c);
uint normal_to_rgb(const vertex &n);
color rgb_to_color(const uint rgb);
vertex rgb_to_normal(const uint n);
color hsv_to_rgb(const int hue, const double saturation, const double value);

/* bounding box related functions */
bounds calc_bounding_box(const std::vector<face *> &faces);
void intersect_bounding_boxes(const bounds &b1, const bounds &b2,
                              bounds &newbounds);
void list_bounds(bounds &newbounds, int nverts, const vertex *vertices...);

/* geometry intersections */
bool ray_AABB_intersect(const bounds &AABB, const ray &r);
face *ray_faces_intersect(const std::vector<face *> &faces, const ray &r,
                          bool lazy, vertex *tuv);
bool ray_sphere_intersect(const ray &r, const double rad,
                          double &t); // intersects a ray with a sphere centered
                                      // on the origin. Assumes ray direction is
                                      // normalized

/* misc geometry functions */
point_2d<double> random_point_in_unit_circle();

#endif /* defined(__Renderer__geom__) */
