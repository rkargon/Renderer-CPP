#include "sampling.h"

#include <cmath>
#include <random>

// // Samples uniformly from a hemisphere around a normal vector
// ray_sample sample_uniform_hemisphere(const Vector3f &n) {
//   Vector3f v = sample_uniform_sphere();
//   if (v.dot(n) < 0) {
//     v = -v;
//   }
//   return {v, 1.0 / (2 * M_PI)};
// }

direction_sample uniform_unit_sphere() {
  double longitude =
      std::uniform_real_distribution<double>(0, 2 * M_PI)(generator);
  double u = std::uniform_real_distribution<double>(-1, 1)(generator);
  double latitude = std::asin(u);

  double cos_long = std::cos(longitude);
  double sin_long = std::sin(longitude);
  double cos_lat = std::cos(latitude);
  double sin_lat = std::sin(latitude);
  return {vertex(cos_lat * cos_long, cos_lat * sin_long, sin_lat),
          1.0 / (4 * M_PI)};
}

direction_sample uniform_unit_hemisphere(const vertex &normal) {
  auto sphere_sample = uniform_unit_sphere();
  if (glm::dot(sphere_sample.first, normal) < 0) {
    sphere_sample.first *= -1;
  }
  sphere_sample.second = 1.0 / (2 * M_PI);
  return sphere_sample;
}
